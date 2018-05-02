WOMTOOL_VERSION=31

if [ ! -f womtool.jar ]
then
  curl -fsSL https://github.com/broadinstitute/cromwell/releases/download/31/womtool-${WOMTOOL_VERSION}.jar -o womtool.jar
fi

cd workflows/
for workflow in */
do
  cd ${workflow}
  for version in */
  do
    cd ${version}
    if [ ! -f ${workflow%/}.wdl ]
    then
      cd ..
      continue
    fi
    echo "Generating inputs.json file for workflow ${workflow%/} version ${version%/}"
    sed 's|https://raw.githubusercontent.com/labbcb/rnnr/master|../../..|g' < ${workflow%/}.wdl > tmp.wdl
    java -jar ../../../womtool.jar inputs tmp.wdl > inputs.json
    rm tmp.wdl
    cd ..
  done
  cd ..
done
