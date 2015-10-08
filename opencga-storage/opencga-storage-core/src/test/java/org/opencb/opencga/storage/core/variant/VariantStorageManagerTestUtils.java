package org.opencb.opencga.storage.core.variant;

import org.junit.*;
import org.junit.rules.ExpectedException;
import org.opencb.biodata.formats.io.FileFormatException;
import org.opencb.commons.test.GenericTest;
import org.opencb.datastore.core.ObjectMap;
import org.opencb.opencga.core.common.IOUtils;
import org.opencb.opencga.storage.core.StorageManagerException;
import org.opencb.opencga.storage.core.StudyConfiguration;
import org.opencb.opencga.storage.core.variant.annotation.VariantAnnotationManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

/**
 * Created by jacobo on 31/05/15.
 */
@Ignore
public abstract class VariantStorageManagerTestUtils extends GenericTest {


    public static final String VCF_TEST_FILE_NAME = "10k.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz";
    public static final String SMALL_VCF_TEST_FILE_NAME = "variant-test-file.vcf.gz";
    public static final String VCF_CORRUPTED_FILE_NAME = "variant-test-file-corrupted.vcf";
    public static final int NUM_VARIANTS = 9792;
    public static final int STUDY_ID = 1;
    public static final String STUDY_NAME = "1000g";
    public static final String DB_NAME = "opencga_variants_test";

    protected static URI inputUri;
    protected static URI smallInputUri;
    protected static URI corruptedInputUri;
    protected static URI outputUri;
    protected VariantStorageManager variantStorageManager;
    public static Logger logger;

    @Rule
    public ExpectedException thrown = ExpectedException.none();

    @BeforeClass
    public static void _beforeClass() throws Exception {
//        System.setProperty(org.slf4j.impl.SimpleLogger.DEFAULT_LOG_LEVEL_KEY, "debug");
        Path rootDir = getTmpRootDir();
        if (rootDir.toFile().exists()) {
            IOUtils.deleteDirectory(rootDir);
            Files.createDirectories(rootDir);
        }
        Path inputPath = rootDir.resolve(VCF_TEST_FILE_NAME);
        Path smallInputPath = rootDir.resolve(SMALL_VCF_TEST_FILE_NAME);
        Path corruptedInputPath = rootDir.resolve(VCF_CORRUPTED_FILE_NAME);
        Files.copy(VariantStorageManagerTest.class.getClassLoader().getResourceAsStream(VCF_TEST_FILE_NAME), inputPath, StandardCopyOption.REPLACE_EXISTING);
        Files.copy(VariantStorageManagerTest.class.getClassLoader().getResourceAsStream(SMALL_VCF_TEST_FILE_NAME), smallInputPath, StandardCopyOption.REPLACE_EXISTING);
        Files.copy(VariantStorageManagerTest.class.getClassLoader().getResourceAsStream(VCF_CORRUPTED_FILE_NAME), corruptedInputPath, StandardCopyOption.REPLACE_EXISTING);

        inputUri = inputPath.toUri();
        smallInputUri = smallInputPath.toUri();
        corruptedInputUri = corruptedInputPath.toUri();
        outputUri = rootDir.toUri();
        logger = LoggerFactory.getLogger(VariantStorageManagerTest.class);

    }

    public static URI getResourceUri(String resourceName) throws IOException {
        Path rootDir = getTmpRootDir();
        Path resourcePath = rootDir.resolve(resourceName);
        Files.copy(VariantStorageManagerTest.class.getClassLoader().getResourceAsStream(resourceName), resourcePath, StandardCopyOption.REPLACE_EXISTING);
        return resourcePath.toUri();
    }

    protected static Path getTmpRootDir() throws IOException {
        Path rootDir = Paths.get("/tmp", "VariantStorageManagerTest");
        Files.createDirectories(rootDir);
        return rootDir;
    }

    @Before
    public void before() throws Exception {
        clearDB(DB_NAME);
    }

    @Before
    public final void _before() throws Exception {
        variantStorageManager = getVariantStorageManager();
    }

    protected abstract VariantStorageManager getVariantStorageManager() throws Exception;
    protected abstract void clearDB(String dbName) throws Exception;


    /* ---------------------------------------------------- */
    /* Static methods to run a simple ETL to index Variants */
    /* ---------------------------------------------------- */

    /**
     * Simple class to store the output URIs generated by the ETL
     */
    public static class ETLResult {

        public URI extractResult;
        public URI preTransformResult;
        public URI transformResult;
        public URI postTransformResult;
        public URI preLoadResult;
        public URI loadResult;
        //        public URI postLoadResult;
    }

    public static ETLResult runETL(VariantStorageManager variantStorageManager, ObjectMap options)
            throws IOException, FileFormatException, StorageManagerException {
        return runETL(variantStorageManager, options, true, true, true);
    }

    public static ETLResult runETL(VariantStorageManager variantStorageManager, ObjectMap options,
                                   boolean doExtract,
                                   boolean doTransform,
                                   boolean doLoad)
            throws IOException, FileFormatException, StorageManagerException {
        return runETL(variantStorageManager, inputUri, outputUri, options, options, options, options, options, options, options, doExtract, doTransform, doLoad);
    }

    public static ETLResult runDefaultETL(VariantStorageManager variantStorageManager, StudyConfiguration studyConfiguration)
            throws URISyntaxException, IOException, FileFormatException, StorageManagerException {
        return runDefaultETL(inputUri, variantStorageManager, studyConfiguration);
    }

    public static ETLResult runDefaultETL(URI inputUri, VariantStorageManager variantStorageManager, StudyConfiguration studyConfiguration)
            throws URISyntaxException, IOException, FileFormatException, StorageManagerException {
        return runDefaultETL(inputUri, variantStorageManager, studyConfiguration, new ObjectMap());
    }

    public static ETLResult runDefaultETL(URI inputUri, VariantStorageManager variantStorageManager,
                                          StudyConfiguration studyConfiguration, ObjectMap params)
            throws URISyntaxException, IOException, FileFormatException, StorageManagerException {

        ObjectMap extractParams = new ObjectMap(params);

        ObjectMap preTransformParams = new ObjectMap(params);
        preTransformParams.put(VariantStorageManager.Options.STUDY_CONFIGURATION.key(), studyConfiguration);
        preTransformParams.putIfAbsent(VariantStorageManager.Options.FILE_ID.key(), 6);

        ObjectMap transformParams = new ObjectMap(params);
        transformParams.put(VariantStorageManager.Options.STUDY_CONFIGURATION.key(), studyConfiguration);
        transformParams.putIfAbsent(VariantStorageManager.Options.INCLUDE_GENOTYPES.key(), true);
        transformParams.putIfAbsent(VariantStorageManager.Options.FILE_ID.key(), 6);

        ObjectMap postTransformParams = new ObjectMap(params);

        ObjectMap preLoadParams = new ObjectMap(params);
        preLoadParams.put(VariantStorageManager.Options.STUDY_CONFIGURATION.key(), studyConfiguration);

        ObjectMap loadParams = new ObjectMap(params);
        loadParams.put(VariantStorageManager.Options.STUDY_CONFIGURATION.key(), studyConfiguration);
        loadParams.putIfAbsent(VariantStorageManager.Options.INCLUDE_GENOTYPES.key(), true);
        loadParams.putIfAbsent(VariantStorageManager.Options.FILE_ID.key(), 6);
        loadParams.putIfAbsent(VariantStorageManager.Options.DB_NAME.key(), DB_NAME);

        ObjectMap postLoadParams = new ObjectMap(params);
        postLoadParams.put(VariantStorageManager.Options.STUDY_CONFIGURATION.key(), studyConfiguration);
        postLoadParams.putIfAbsent(VariantStorageManager.Options.DB_NAME.key(), DB_NAME);
        postLoadParams.putIfAbsent(VariantStorageManager.Options.FILE_ID.key(), 6);
        postLoadParams.putIfAbsent(VariantStorageManager.Options.ANNOTATE.key(), true);
        postLoadParams.putIfAbsent(VariantAnnotationManager.SPECIES, "hsapiens");
        postLoadParams.putIfAbsent(VariantAnnotationManager.ASSEMBLY, "GRc37");
        postLoadParams.putIfAbsent(VariantStorageManager.Options.CALCULATE_STATS.key(), true);

        return runETL(variantStorageManager, inputUri, outputUri, extractParams, preTransformParams, transformParams, postTransformParams, preLoadParams, loadParams, postLoadParams, true, true, true);
    }

    public static ETLResult runETL(VariantStorageManager variantStorageManager, URI inputUri, URI outputUri,
                                   ObjectMap extractParams,
                                   ObjectMap preTransformParams, ObjectMap transformParams, ObjectMap postTransformParams,
                                   ObjectMap preLoadParams, ObjectMap loadParams, ObjectMap postLoadParams,
                                   boolean doExtract,
                                   boolean doTransform,
                                   boolean doLoad)
            throws IOException, FileFormatException, StorageManagerException {
        ETLResult etlResult = new ETLResult();

        if (doExtract) {
            variantStorageManager.getConfiguration().getStorageEngine(variantStorageManager.getStorageEngineId()).getVariant().setOptions(extractParams);
            inputUri = variantStorageManager.extract(inputUri, outputUri);
            etlResult.extractResult = inputUri;
        }

        if (doTransform) {
            variantStorageManager.getConfiguration().getStorageEngine(variantStorageManager.getStorageEngineId()).getVariant().setOptions(preTransformParams);
            inputUri = variantStorageManager.preTransform(inputUri);
            etlResult.preTransformResult = inputUri;
            Assert.assertTrue("Intermediary file " + inputUri + " does not exist", Paths.get(inputUri).toFile().exists());


            variantStorageManager.getConfiguration().getStorageEngine(variantStorageManager.getStorageEngineId()).getVariant().setOptions(transformParams);
            inputUri = variantStorageManager.transform(inputUri, null, outputUri);
            etlResult.transformResult = inputUri;
            Assert.assertTrue("Intermediary file " + inputUri + " does not exist", Paths.get(inputUri).toFile().exists());

            variantStorageManager.getConfiguration().getStorageEngine(variantStorageManager.getStorageEngineId()).getVariant().setOptions(postTransformParams);
            inputUri = variantStorageManager.postTransform(inputUri);
            etlResult.postTransformResult = inputUri;
            Assert.assertTrue("Intermediary file " + inputUri + " does not exist", Paths.get(inputUri).toFile().exists());
        }

        if (doLoad) {

            variantStorageManager.getConfiguration().getStorageEngine(variantStorageManager.getStorageEngineId()).getVariant().setOptions(preLoadParams);
            inputUri = variantStorageManager.preLoad(inputUri, outputUri);
            etlResult.preLoadResult = inputUri;
            Assert.assertTrue("Intermediary file " + inputUri + " does not exist", Paths.get(inputUri).toFile().exists());


            variantStorageManager.getConfiguration().getStorageEngine(variantStorageManager.getStorageEngineId()).getVariant().setOptions(loadParams);
            inputUri = variantStorageManager.load(inputUri);
            etlResult.loadResult = inputUri;
            Assert.assertTrue("Intermediary file " + inputUri + " does not exist", Paths.get(inputUri).toFile().exists());


            variantStorageManager.getConfiguration().getStorageEngine(variantStorageManager.getStorageEngineId()).getVariant().setOptions(postLoadParams);
            variantStorageManager.postLoad(inputUri, outputUri);
        }
        return etlResult;
    }

    protected static StudyConfiguration newStudyConfiguration() {
        return new StudyConfiguration(STUDY_ID, STUDY_NAME);
    }


}
