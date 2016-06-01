package org.opencb.opencga.analysis.execution.plugins.diagnosis;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.commons.lang3.RandomUtils;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.commons.datastore.core.ObjectMap;
import org.opencb.commons.datastore.core.Query;
import org.opencb.commons.datastore.core.QueryOptions;
import org.opencb.opencga.analysis.execution.plugins.OpenCGAAnalysis;
import org.opencb.opencga.catalog.models.DiseasePanel;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBAdaptor;
import org.opencb.opencga.storage.core.variant.adaptors.VariantDBIterator;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * Created by mmedina on 5/31/16.
 */
public class DiagnosisAnalysis extends OpenCGAAnalysis {

    @Override
    public String getIdentifier() {
        return "diagnosis";
    }

    @Override
    public int run() throws Exception {
        long studyId = getStudyId();
        String sessionId = getSessionId();
        List<String> samples = getConfiguration().getAsStringList("samples");
        String outdir = getConfiguration().getString("outdir");
        String panelName = getConfiguration().getString("panel");

        DiseasePanel panel = getCatalogManager().getDiseasePanel(panelName, null, sessionId).first();

        VariantDBAdaptor dbAdaptor = getVariantDBAdaptor(studyId);

        Query query = new Query(VariantDBAdaptor.VariantQueryParams.STUDIES.key(), studyId)
                .append(VariantDBAdaptor.VariantQueryParams.RETURNED_SAMPLES.key(), samples);
        VariantDBIterator iterator = dbAdaptor.iterator(query, new QueryOptions());

        Map<String, List<String>> diagnosticVariants = initMap(samples);
        Map<String, List<String>> secondaryFinding = initMap(samples);
        while (iterator.hasNext()) {
            Variant variant = iterator.next();
            /**
             * DO MAGIC
             */
            if (RandomUtils.nextInt(0, 1000) == 0) {
                diagnosticVariants.get(samples.get(0)).add(variant.toString());
            }

        }

        for (String sample : samples) {
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
