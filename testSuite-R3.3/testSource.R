
  message("STARTING JunctionSeq Post-Build Test (",TEST.ID,") (",date(),")");

  sc <- file(paste0("logs/",TEST.ID,".log"), open = 'w');
  sink(file = sc , type = "message");

  message("STARTING JunctionSeq Post-Build Test (",TEST.ID,") (",date(),")");

  errorFlagFile=paste0("status/",TEST.ID,"-BAD_ERRORS.txt");
  warnFlagFile=paste0("status/",TEST.ID,"-BAD_WARNINGS.txt");
  okFlagFile=paste0("status/",TEST.ID,"-OK_COMPLETED_NO_WARNINGS.txt");

  if(file.exists(errorFlagFile)){ file.remove(errorFlagFile) }
  if(file.exists(warnFlagFile)){ file.remove(warnFlagFile) }
  if(file.exists(okFlagFile)){ file.remove(okFlagFile) }

  no.errors <- TRUE;

  tryCatch({
      source("testFunctions.R");
      source(paste0("",TEST.ID,".R"));
  }, error = function(e){
    no.errors <- FALSE;
    write.table(c(),file=errorFlagFile, quote=F);
    stop(e);
  });

  no.warnings <- is.null(warnings());

  if(no.warnings & no.errors){
    write.table(c(),file=okFlagFile, quote=F);
  } else {
    write.table(c(),file=warnFlagFile, quote=F);
  }
  
  warnings();
  
  message("FINISHED JunctionSeq Post-Build Tests... (",date(),")");

  sink(file = NULL, type = "message");
  close(sc);

  message("FINISHED JunctionSeq Post-Build Tests... (",date(),")");

