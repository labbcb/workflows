#!/usr/bin/env bash
WOMTOOL_VERSION="48"

if [ ! -f womtool.jar ]
then
  echo "Downloading WOMTOOL version ${WOMTOOL_VERSION}"
  curl -fsSL "https://github.com/broadinstitute/cromwell/releases/download/${WOMTOOL_VERSION}/womtool-${WOMTOOL_VERSION}.jar" -o womtool.jar
  if [ $? -ne 0 ]; then
    exit 1
  fi
fi

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
    java -jar ../../womtool.jar inputs ${workflow%/}.wdl > inputs.json
    cd ..
  done
  cd ..
done
