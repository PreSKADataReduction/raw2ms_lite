# Installing `singularity`

Following the [installation tutorial](https://singularity-tutorial.github.io/01-installation/)

```
sudo apt install -y build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools libseccomp-dev wget pkg-config git cryptsetup debootstrap 

sudo apt install -y containernetworking-plugins golang-go

wget http://http.us.debian.org/debian/pool/main/s/singularity-container/singularity-container_3.5.2+ds1-1_amd64.deb

sudo dpkg -i singularity-container_3.5.2+ds1-1_amd64.deb  
```
