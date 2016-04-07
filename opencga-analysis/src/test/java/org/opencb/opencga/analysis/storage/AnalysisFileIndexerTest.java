package org.opencb.opencga.analysis.storage;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;
import org.opencb.biodata.models.variant.VariantSource;
import org.opencb.commons.datastore.core.Query;
import org.opencb.commons.datastore.core.QueryOptions;
import org.opencb.datastore.mongodb.MongoDataStore;
import org.opencb.datastore.mongodb.MongoDataStoreManager;
import org.opencb.opencga.analysis.files.FileMetadataReader;
import org.opencb.opencga.catalog.CatalogManager;
import org.opencb.opencga.catalog.CatalogManagerTest;
import org.opencb.opencga.catalog.config.CatalogConfiguration;
import org.opencb.opencga.catalog.config.Policies;
import org.opencb.opencga.catalog.db.api.CatalogCohortDBAdaptor;
import org.opencb.opencga.catalog.db.api.CatalogFileDBAdaptor;
import org.opencb.opencga.catalog.exceptions.CatalogException;
import org.opencb.opencga.catalog.models.*;
import org.opencb.opencga.catalog.utils.CatalogFileUtils;
import org.opencb.opencga.storage.core.variant.VariantStorageManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.net.URI;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.opencb.biodata.models.variant.StudyEntry.DEFAULT_COHORT;
import static org.opencb.opencga.analysis.storage.AnalysisStorageTestUtil.runStorageJob;
import static org.opencb.opencga.analysis.storage.variant.VariantStorageTest.checkCalculatedAggregatedStats;
import static org.opencb.opencga.analysis.storage.variant.VariantStorageTest.checkCalculatedStats;
import static org.opencb.opencga.storage.core.variant.VariantStorageManagerTestUtils.DB_NAME;
import static org.opencb.opencga.storage.core.variant.VariantStorageManagerTestUtils.getResourceUri;

/**
 * Created by hpccoll1 on 13/07/15.
 */
public class AnalysisFileIndexerTest {

    @Rule
    public ExpectedException thrown = ExpectedException.none();

    private CatalogManager catalogManager;
    private String sessionId;
    private long projectId;
    private long studyId;
    private FileMetadataReader fileMetadataReader;
    private CatalogFileUtils catalogFileUtils;
    private long outputId;
    Logger logger = LoggerFactory.getLogger(AnalysisFileIndexerTest.class);
    private CatalogConfiguration catalogConfiguration;
    private final String userId = "user";
    private final String dbName = DB_NAME;
    private List<File> files = new ArrayList<>();

    public void beforeIndex() throws Exception {
        Path openCGA = AnalysisStorageTestUtil.isolateOpenCGA();
        catalogConfiguration = CatalogConfiguration.load(getClass().getResource("/catalog-configuration.yml").openStream());
//        Properties properties = Config.getCatalogProperties();

        CatalogManagerTest.clearCatalog(catalogConfiguration);
        clearDB(dbName);

        catalogManager = new CatalogManager(catalogConfiguration);
        fileMetadataReader = FileMetadataReader.get(catalogManager);
        catalogFileUtils = new CatalogFileUtils(catalogManager);

        User user = catalogManager.createUser(userId, "User", "user@email.org", "user", "ACME", null).first();
        sessionId = catalogManager.login(userId, "user", "localhost").first().getString("sessionId");
        projectId = catalogManager.createProject(userId, "p1", "p1", "Project 1", "ACME", null, sessionId).first().getId();
        studyId = catalogManager.createStudy(projectId, "s1", "s1", Study.Type.CASE_CONTROL, null, null, "Study 1", null, null, null, null, Collections.singletonMap(File.Bioformat.VARIANT, new DataStore("mongodb", dbName)), null, null, null, sessionId).first().getId();
        outputId = catalogManager.createFolder(studyId, Paths.get("data", "index"), false, null, sessionId).first().getId();
        files.add(create("1-500.filtered.10k.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"));
        files.add(create("501-1000.filtered.10k.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"));
        files.add(create("1001-1500.filtered.10k.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"));
        files.add(create("1501-2000.filtered.10k.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"));
        files.add(create("2001-2504.filtered.10k.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"));

    }

    private void clearDB(String dbName) {
        logger.info("Cleaning MongoDB {}" , dbName);
        MongoDataStoreManager mongoManager = new MongoDataStoreManager("localhost", 27017);
        MongoDataStore mongoDataStore = mongoManager.get(dbName);
        mongoManager.drop(dbName);
    }

    public File create(String resourceName) throws IOException, CatalogException {
        File file;
        URI uri = getResourceUri(resourceName);
        file = fileMetadataReader.create(studyId, uri, "data/vcfs/", "", true, null, sessionId).first();
        catalogFileUtils.upload(uri, file, null, sessionId, false, false, true, false, Long.MAX_VALUE);
        return catalogManager.getFile(file.getId(), sessionId).first();
    }

    @Test
    public void testIndexWithStats() throws Exception {
        beforeIndex();

        AnalysisFileIndexer analysisFileIndexer = new AnalysisFileIndexer(catalogManager);
        QueryOptions queryOptions = new QueryOptions(VariantStorageManager.Options.ANNOTATE.key(), false);
        runStorageJob(catalogManager, analysisFileIndexer.index((int) files.get(0).getId(), (int) outputId, sessionId, queryOptions).first(), logger, sessionId);
        assertEquals(500, getDefaultCohort().getSamples().size());
        assertEquals(Cohort.CohortStatus.NONE, getDefaultCohort().getStatus().getStatus());
        assertNotNull(catalogManager.getFile(files.get(0).getId(), sessionId).first().getStats().get(FileMetadataReader.VARIANT_STATS));

        runStorageJob(catalogManager, analysisFileIndexer.index((int) files.get(1).getId(), (int) outputId, sessionId, queryOptions).first(), logger, sessionId);
        assertEquals(1000, getDefaultCohort().getSamples().size());
        assertEquals(Cohort.CohortStatus.NONE, getDefaultCohort().getStatus().getStatus());
        assertNotNull(catalogManager.getFile(files.get(1).getId(), sessionId).first().getStats().get(FileMetadataReader.VARIANT_STATS));

        queryOptions.put(VariantStorageManager.Options.CALCULATE_STATS.key(), true);
        runStorageJob(catalogManager, analysisFileIndexer.index((int) files.get(2).getId(), (int) outputId, sessionId, queryOptions).first(), logger, sessionId);
        assertEquals(1500, getDefaultCohort().getSamples().size());
        assertEquals(Cohort.CohortStatus.READY, getDefaultCohort().getStatus().getStatus());
        checkCalculatedStats(Collections.singletonMap(DEFAULT_COHORT, catalogManager.getAllCohorts(studyId,
                new Query(CatalogCohortDBAdaptor.QueryParams.NAME.key(), DEFAULT_COHORT), new QueryOptions(), sessionId).first()),
                catalogManager, dbName, sessionId);
        assertNotNull(catalogManager.getFile(files.get(2).getId(), sessionId).first().getStats().get(FileMetadataReader.VARIANT_STATS));

        queryOptions.put(VariantStorageManager.Options.CALCULATE_STATS.key(), false);
        runStorageJob(catalogManager, analysisFileIndexer.index((int) files.get(3).getId(), (int) outputId, sessionId, queryOptions).first(), logger, sessionId);
        assertEquals(2000, getDefaultCohort().getSamples().size());
        assertEquals(Cohort.CohortStatus.INVALID, getDefaultCohort().getStatus().getStatus());
        assertNotNull(catalogManager.getFile(files.get(3).getId(), sessionId).first().getStats().get(FileMetadataReader.VARIANT_STATS));

        queryOptions.put(VariantStorageManager.Options.CALCULATE_STATS.key(), true);
        runStorageJob(catalogManager, analysisFileIndexer.index((int) files.get(4).getId(), (int) outputId, sessionId, queryOptions).first(), logger, sessionId);
        assertEquals(2504, getDefaultCohort().getSamples().size());
        assertEquals(Cohort.CohortStatus.READY, getDefaultCohort().getStatus().getStatus());
        assertNotNull(catalogManager.getFile(files.get(4).getId(), sessionId).first().getStats().get(FileMetadataReader.VARIANT_STATS));

        checkCalculatedStats(Collections.singletonMap(DEFAULT_COHORT, catalogManager.getAllCohorts(studyId,
                new Query(CatalogCohortDBAdaptor.QueryParams.NAME.key(), DEFAULT_COHORT), new QueryOptions(), sessionId).first()),
                catalogManager, dbName, sessionId);
    }
    

    private Cohort getDefaultCohort() throws CatalogException {
        return catalogManager.getAllCohorts(studyId, new Query(CatalogCohortDBAdaptor.QueryParams.NAME.key(), DEFAULT_COHORT),
                new QueryOptions(), sessionId).first();
    }

    @Test
    public void testIndexBySteps() throws Exception {
        beforeIndex();
        AnalysisFileIndexer analysisFileIndexer = new AnalysisFileIndexer(catalogManager);
        QueryOptions queryOptions = new QueryOptions(VariantStorageManager.Options.ANNOTATE.key(), false).append(VariantStorageManager.Options.CALCULATE_STATS.key(), true);

        queryOptions.append(AnalysisFileIndexer.TRANSFORM, true);
        queryOptions.append(AnalysisFileIndexer.LOAD, false);

        //Create transform index job
        Job job = analysisFileIndexer.index((int) files.get(0).getId(), (int) outputId, sessionId, queryOptions).first();
        assertEquals(Index.IndexStatus.TRANSFORMING, catalogManager.getFile(files.get(0).getId(), sessionId).first().getIndex().getStatus().getStatus());

        //Run transform index job
        job = runStorageJob(catalogManager, job, logger, sessionId);
        assertEquals(500, getDefaultCohort().getSamples().size());
        assertEquals(Cohort.CohortStatus.NONE, getDefaultCohort().getStatus().getStatus());
        assertEquals(Index.IndexStatus.TRANSFORMED, catalogManager.getFile(files.get(0).getId(), sessionId).first().getIndex().getStatus().getStatus());
        assertEquals(job.getAttributes().get(Job.TYPE), Job.Type.INDEX.toString());

        //Get transformed file
        Query searchQuery = new Query(CatalogFileDBAdaptor.QueryParams.ID.key(), job.getOutput())
                .append(CatalogFileDBAdaptor.QueryParams.NAME.key(), "~variants.(json|avro)");
        File transformedFile = catalogManager.getAllFiles(studyId, searchQuery, new QueryOptions(), sessionId).first();
        assertEquals(job.getId(), transformedFile.getJobId());
        File file = catalogManager.getFile(files.get(0).getId(), sessionId).first();
        assertNotNull(file.getStats().get(FileMetadataReader.VARIANT_STATS));

        //Create load index job
        queryOptions.append(AnalysisFileIndexer.TRANSFORM, false);
        queryOptions.append(AnalysisFileIndexer.LOAD, true);
        job = analysisFileIndexer.index((int) transformedFile.getId(), (int) outputId, sessionId, queryOptions).first();
        assertEquals(Index.IndexStatus.LOADING, catalogManager.getFile(files.get(0).getId(), sessionId).first().getIndex().getStatus().getStatus());
        assertEquals(job.getAttributes().get(Job.TYPE), Job.Type.INDEX.toString());

        //Run load index job
        runStorageJob(catalogManager, job, logger, sessionId);
        assertEquals(500, getDefaultCohort().getSamples().size());
        assertEquals(Cohort.CohortStatus.READY, getDefaultCohort().getStatus().getStatus());
        assertEquals(Index.IndexStatus.READY, catalogManager.getFile(files.get(0).getId(), sessionId).first().getIndex().getStatus().getStatus());
        Cohort defaultCohort = catalogManager.getAllCohorts(studyId,
                new Query(CatalogCohortDBAdaptor.QueryParams.NAME.key(), DEFAULT_COHORT), new QueryOptions(), sessionId).first();
        checkCalculatedStats(Collections.singletonMap(DEFAULT_COHORT, defaultCohort), catalogManager, dbName, sessionId);
    }
    
    public void beforeAggregatedIndex(String file, VariantSource.Aggregation aggregation) throws Exception {
        Path openCGA = AnalysisStorageTestUtil.isolateOpenCGA();
        catalogConfiguration = CatalogConfiguration.load(getClass().getResource("/catalog-configuration.yml").openStream());
        //Properties properties = Config.getCatalogProperties();

        CatalogManagerTest.clearCatalog(catalogConfiguration);
        clearDB(dbName);

        catalogManager = new CatalogManager(catalogConfiguration);
        fileMetadataReader = FileMetadataReader.get(catalogManager);
        catalogFileUtils = new CatalogFileUtils(catalogManager);
        Policies policies = new Policies();
        policies.setUserCreation(Policies.UserCreation.ALWAYS);

        User user = catalogManager.createUser(userId, "User", "user@email.org", "user", "ACME", null).first();
        sessionId = catalogManager.login(userId, "user", "localhost").first().getString("sessionId");
        projectId = catalogManager.createProject(userId, "p1", "p1", "Project 1", "ACME", null, sessionId).first().getId();
        studyId = catalogManager.createStudy(projectId, "s1", "s1", Study.Type.CASE_CONTROL, null, null, "Study 1", null,
                null, null, null, Collections.singletonMap(File.Bioformat.VARIANT, new DataStore("mongodb", dbName)), null,
                Collections.singletonMap(VariantStorageManager.Options.AGGREGATED_TYPE.key(), VariantSource.Aggregation.BASIC), 
                null, sessionId).first().getId();
        outputId = catalogManager.createFolder(studyId, Paths.get("data", "index"), false, null, sessionId).first().getId();
        files.add(create("variant-test-aggregated-file.vcf.gz"));
    }
    
    @Test
    public void testIndexWithAggregatedStats() throws Exception {
        beforeAggregatedIndex("variant-test-aggregated-file.vcf.gz", VariantSource.Aggregation.BASIC);

        AnalysisFileIndexer analysisFileIndexer = new AnalysisFileIndexer(catalogManager);
        QueryOptions queryOptions = new QueryOptions(VariantStorageManager.Options.ANNOTATE.key(), false)
                .append(VariantStorageManager.Options.AGGREGATED_TYPE.key(), VariantSource.Aggregation.BASIC);

        queryOptions.put(VariantStorageManager.Options.CALCULATE_STATS.key(), true);
        runStorageJob(catalogManager, analysisFileIndexer.index((int) files.get(0).getId(), (int) outputId, sessionId, queryOptions).first(), logger, sessionId);
        assertEquals(0, getDefaultCohort().getSamples().size());
        assertEquals(Cohort.CohortStatus.READY, getDefaultCohort().getStatus().getStatus());
        checkCalculatedAggregatedStats(Collections.singletonMap(DEFAULT_COHORT, catalogManager.getAllCohorts(studyId,
                new Query(CatalogCohortDBAdaptor.QueryParams.NAME.key(), DEFAULT_COHORT), new QueryOptions(), sessionId).first()), dbName
        );
    }
    
}


