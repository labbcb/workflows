#!/usr/bin/env bash
cd tools/
for tool in */
do
  cd ${tool}
  for version in */
  do
    cd ${version}
    if [ ! -f Dockerfile ]
    then
      cd ..
      continue
    fi
    docker build --tag "welliton/${tool%/}:${version%/}" .
    if [ $? -ne 0 ]; then
        exit 1
    fi
    docker push "welliton/${tool%/}:${version%/}"
    if [ $? -ne 0 ]; then
        exit 1
    fi
    cd ..
  done
  cd ..
done
