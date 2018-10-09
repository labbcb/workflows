#!/usr/bin/env bash
WOMTOOL_VERSION="35"

if [ ! -f womtool.jar ]
then
  echo "Download WOMTOOL version ${WOMTOOL_VERSION}"
  curl -fsSL "https://github.com/broadinstitute/cromwell/releases/download/${WOMTOOL_VERSION}/womtool-${WOMTOOL_VERSION}.jar" -o womtool.jar
  if [ $? -ne 0 ]; then
    exit 1
  fi
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
    sed 's|https://raw.githubusercontent.com/labbcb/workflows/master|../../..|g' < ${workflow%/}.wdl > tmp.wdl
    java -jar ../../../womtool.jar inputs tmp.wdl > inputs.json
    rm tmp.wdl
    cd ..
  done
  cd ..
done
