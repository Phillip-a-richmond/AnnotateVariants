BIN_VERSION="1.1.0"

sudo apt -y update
sudo apt-get -y install docker.io
sudo docker pull google/deepvariant:"${BIN_VERSION}"
