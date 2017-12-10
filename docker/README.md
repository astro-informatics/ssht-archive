# Docker container for pyssht

To build the Docker image:
```
docker build -t='eiffl/pyssht:2.0' .
```
To push the image to Dockerhub
```
docker push 'eiffl/pyssht:2.0'
```
To use pyssht within the Docker container:
```
docker run -it "eiffl/pyssht" /usr/bin/python
```
