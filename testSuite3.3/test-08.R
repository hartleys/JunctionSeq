
runJS();

buildAllPlots(jscs, outfile.prefix = outputDir("std"), plotting.device.params = list(res = 72));
buildAllPlots(jscs, outfile.prefix = outputDir("GL3"), gene.list = c(g1,g2,g3), plotting.device.params = list(res = 72));
buildAllPlots(jscs, outfile.prefix = outputDir("G1"), gene.list = c(g1), plotting.device.params = list(res = 72));

buildAllPlots(jscs, outfile.prefix = outputDir("std-LongName"), minimalImageFilenames = FALSE, plotting.device.params = list(res = 72));
buildAllPlots(jscs, outfile.prefix = outputDir("GL3-LongName"), gene.list = c(g1,g2,g3), minimalImageFilenames = FALSE, plotting.device.params = list(res = 72));
buildAllPlots(jscs, outfile.prefix = outputDir("G1-LongName"), gene.list = c(g1), minimalImageFilenames = FALSE, plotting.device.params = list(res = 72));



