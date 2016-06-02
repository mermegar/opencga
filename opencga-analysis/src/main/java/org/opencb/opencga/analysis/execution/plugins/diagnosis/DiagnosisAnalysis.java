package org.opencb.opencga.analysis.execution.plugins.diagnosis;

import org.apache.commons.lang3.RandomUtils;
import org.apache.commons.lang3.StringUtils;
import org.opencb.biodata.models.core.Gene;
import org.opencb.biodata.models.core.Region;
import org.opencb.biodata.models.variant.StudyEntry;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.biodata.models.variant.avro.ConsequenceType;
import org.opencb.cellbase.core.client.CellBaseClient;
import org.opencb.commons.datastore.core.*;
import org.opencb.opencga.analysis.execution.plugins.OpenCGAAnalysis;
import org.opencb.opencga.catalog.db.api.CatalogSampleDBAdaptor;
import org.opencb.opencga.catalog.exceptions.CatalogException;
import org.opencb.opencga.catalog.models.DiseasePanel;
import org.opencb.opencga.catalog.models.Sample;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.nio.file.Paths;
import java.util.*;

/**
 * Created by mmedina on 5/31/16.
 */
public class DiagnosisAnalysis extends OpenCGAAnalysis {

    public static final String IDENTIFIER = "diagnosis";
    public static final String PARAM_SAMPLES = "samples";
    public static final String PARAM_PANEL = "panel";
    public static final Set<String> REFERENCE_GENOTYPES = new HashSet<>();

    static {
        REFERENCE_GENOTYPES.add("0/0");
        REFERENCE_GENOTYPES.add("0|0");
        REFERENCE_GENOTYPES.add("0");
        REFERENCE_GENOTYPES.add(".");
        REFERENCE_GENOTYPES.add("./.");
    }

    @Override
    public String getIdentifier() {
        return IDENTIFIER;
    }

    @Override
    public int run() throws Exception {
        long studyId = getStudyId();
        String sessionId = getSessionId();
        List<String> samples = getConfiguration().getAsStringList("samples");
        String outdir = getConfiguration().getString("outdir");
        String panelName = getConfiguration().getString("panel");

        // Get Panel
        DiseasePanel panel = getCatalogManager().getDiseasePanel(panelName, null, sessionId).first();
        Set<String> panelVariantsSet = new HashSet<>(panel.getVariants());
        List<Region> panelRegionsList = new ArrayList<>(panel.getGenes().size() + panel.getRegions().size());

        // Get regions
        for (String region : panel.getRegions()) {
            panelRegionsList.add(Region.parseRegion(region));
        }

        CellBaseClient cellbase = getCellBaseClient("hsapiens");
        QueryResponse<Gene> geneQueryResponse = cellbase.getInfo(CellBaseClient.Category.feature,
                CellBaseClient.SubCategory.gene,
                String.join(",", panel.getGenes()),
                new QueryOptions("include", "chromosome,start,end"));

        System.out.println(cellbase.getLastQuery());
        for (QueryResult<Gene> queryResult : geneQueryResponse.getResponse()) {
            for (Gene gene : queryResult.getResult()) {
                System.out.println(gene);
                panelRegionsList.add(new Region(gene.getChromosome(), gene.getStart(), gene.getEnd()));
            }
        }

        // Get SampleNames
        List<Long> sampleIds = new ArrayList<>(samples.size());
        List<String> sampleNames = new ArrayList<>(samples.size());
        for (String sampleId : samples) {
            //TODO: Search by Name if not a number missing
            Sample sample = getCatalogManager().getSample(Long.parseLong(sampleId), null, sessionId).first();
            sampleIds.add(sample.getId());
            sampleNames.add(sample.getName());
        }


        // Build query
        Query query = new Query(VariantDBAdaptor.VariantQueryParams.STUDIES.key(), studyId)
                .append(VariantDBAdaptor.VariantQueryParams.RETURNED_SAMPLES.key(), sampleIds)
                .append(VariantDBAdaptor.VariantQueryParams.RETURNED_STUDIES.key(), studyId);
        VariantDBAdaptor dbAdaptor = getVariantDBAdaptor(studyId);
        VariantDBIterator iterator = dbAdaptor.iterator(query, new QueryOptions());

        Map<String, List<String>> diagnosticVariants = initMap(sampleNames);
        Map<String, List<String>> secondaryFinding = initMap(sampleNames);

        while (iterator.hasNext()) {
            Variant variant = iterator.next();
            String chr = variant.getChromosome();
            Integer start = variant.getStart();
            Integer end = variant.getEnd();


            // A variant will be candidate to be diagnostic if is one of the variants in the panel
            boolean diagnosticCandidate = panelVariantsSet.contains(variant.toString());

            Set<String> genesInVariant = new HashSet<>();
            for (ConsequenceType consequenceType : variant.getAnnotation().getConsequenceTypes()) {
                genesInVariant.add(consequenceType.getGeneName());
                genesInVariant.add(consequenceType.getEnsemblGeneId());
            }
            // A variant will be candidate to be secondary finding if the variant has TraitAssociation
            // AND belongs to one gene or region of the panel
            boolean secondaryCandidate = false;
            if (variant.getAnnotation().getVariantTraitAssociation() != null) {
                for (Region region : panelRegionsList) {
                    if (region.overlaps(chr, start, end)) {
                        secondaryCandidate = true;
                        break;
                    }
                }
            }

            if (diagnosticCandidate || secondaryCandidate) {
                StudyEntry studyEntry = variant.getStudies().get(0);
                for (String sampleName : sampleNames) {
                    String gt = studyEntry.getSampleData(sampleName, "GT");
                    if (!REFERENCE_GENOTYPES.contains(gt)) {
                        if (diagnosticCandidate) {
                            diagnosticVariants.get(sampleName).add(variant.toString());
                        } else if (secondaryCandidate) {
                            secondaryFinding.get(sampleName).add(variant.toString());
                        }
                    }
                }
            }
        }

        // Print results
        for (String sample : sampleNames) {
            try (OutputStream os = new FileOutputStream(Paths.get(outdir).resolve(getOutputFileName(sample)).toFile())) {
                ObjectMap result = new ObjectMap();
                List<String> diagnostic = diagnosticVariants.get(sample);
                if (diagnostic != null) {
                    result.put("diagnostic", diagnostic);
                } else {
                    result.put("diagnostic", Collections.emptyList());
                }

                List<String> secondary = secondaryFinding.get(sample);
                if (secondary != null) {
                    result.put("secondary", secondary);
                } else {
                    result.put("secondary", Collections.emptyList());
                }

                os.write(result.toJson().getBytes());
                os.write('\n');
            }
        }

        return 0;
    }

    private HashMap<String, List<String>> initMap(List<String> samples) {
        HashMap<String, List<String>> map = new HashMap<>();
        for (String sample : samples) {
            map.put(sample, new LinkedList<>());
        }
        return map;
    }

    private String getOutputFileName(String sample) {
        return "diagnosis_" + sample + ".json";
    }
}
