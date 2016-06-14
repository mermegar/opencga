/*
 * Copyright 2015 OpenCB
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.opencb.opencga.app.cli.admin;


import org.opencb.opencga.catalog.exceptions.CatalogException;
import org.opencb.opencga.server.rest.RestServer;

import javax.ws.rs.client.Client;
import javax.ws.rs.client.ClientBuilder;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.Response;
import java.nio.file.Paths;

/**
 * Created by imedina on 02/03/15.
 */
public class ServerCommandExecutor extends AdminCommandExecutor {

    private AdminCliOptionsParser.ServerCommandOptions serverCommandOptions;

    public ServerCommandExecutor(AdminCliOptionsParser.ServerCommandOptions serverCommandOptions) {
        super(serverCommandOptions.commonOptions);
        this.serverCommandOptions = serverCommandOptions;
    }

    @Override
    public void execute() throws Exception {
        logger.debug("Executing variant command line");

        String subCommandString = serverCommandOptions.getParsedSubCommand();
        switch (subCommandString) {
            case "rest":
                rest();
                break;
            case "grpc":
                grpc();
                break;
            default:
                logger.error("Subcommand not valid");
                break;
        }

    }

    private void rest() throws Exception {
        if (serverCommandOptions.restServerCommandOptions.start) {
//            StorageConfiguration storageConfiguration = configuration;
//            if (StringUtils.isNotEmpty(restCommandOptions.restStartCommandOptions.commonOptions.conf)) {
//                Path path = Paths.get(restCommandOptions.restStartCommandOptions.commonOptions.conf);
//                if (Files.exists(path)) {
//                    storageConfiguration = StorageConfiguration.load(Files.newInputStream(path));
//                }
//            }

//            if (StringUtils.isNotEmpty(restCommandOptions.restStartCommandOptions.commonOptions.storageEngine)) {
//                storageConfiguration.setDefaultStorageEngineId(restCommandOptions.restStartCommandOptions.commonOptions.storageEngine);
//            }

            // Server crated and started
            RestServer server = new RestServer(Paths.get(appHome).resolve("conf"));
            server.start();
            server.blockUntilShutdown();
            logger.info("Shutting down OpenCGA Storage REST server");
        }

        if (serverCommandOptions.restServerCommandOptions.stop) {
//            if (serverCommandOptions.restStopCommandOptions.port > 0) {
//                port = restCommandOptions.restStopCommandOptions.port;
//            }

            Client client = ClientBuilder.newClient();
            WebTarget target = client.target("http://localhost:" + 9090)
                    .path("opencga")
                    .path("webservices")
                    .path("rest")
                    .path("admin")
                    .path("stop");
            Response response = target.request().get();
            logger.info(response.toString());
        }
    }

    private void grpc() throws CatalogException {

    }

}
