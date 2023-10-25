# Workflow Description Language (WDL) for rmats2sashimiplot

* [rmats2sashimiplot.wdl](rmats2sashimiplot.wdl) is a WDL ([https://openwdl.org](https://openwdl.org)) implementation of rmats2sashimiplot

## Run on AnVIL

* [https://anvil.terra.bio](https://anvil.terra.bio)
* Create a new workflow using the AnVIL website and paste or upload [rmats2sashimiplot.wdl](rmats2sashimiplot.wdl)

## Run with cromwell

* [https://github.com/broadinstitute/cromwell](https://github.com/broadinstitute/cromwell)
* Create an input config file
  + `java -jar "${WOMTOOL_JAR}" inputs --optional-inputs true rmats2sashimiplot.wdl > rmats2sashimiplot_inputs.json`
  + Edit `rmats2sashimiplot_inputs.json`
* Validate
  + `java -jar "${WOMTOOL_JAR}" validate --inputs rmats2sashimiplot_inputs.json rmats2sashimiplot.wdl`
* Run
  + `java -jar "${CROMWELL_JAR}" run --inputs rmats2sashimiplot_inputs.json rmats2sashimiplot.wdl`
