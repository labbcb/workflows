cd tools/
for tool in */
do
  cd ${tool}
  for version in */
  do
    cd ${version}
    if [ ! -f Dockerfile ]
    then
      continue
    fi
    docker build --tag "welliton/${tool%/}:${version%/}" .
    docker push "welliton/${tool%/}:${version%/}"
    cd ..
  done
done
