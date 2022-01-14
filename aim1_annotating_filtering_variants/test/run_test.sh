# Export the cromwell server path in the below line
# export CROMWELL=Enter the CROMWELL Installation Path:$CROMWELL
# Configure the cromwell.conf file
java -Dconfig.file=cromwell.conf -jar $CROMWELL run -i input.json ../wdl/Annotate_and_Filter_Variants.wdl
