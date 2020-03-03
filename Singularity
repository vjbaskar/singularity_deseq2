BootStrap: docker
From: ubuntu:18.04


%files
	install_packages.R /tmp
	R-3.5.0.tar.gz /tmp
	deseq2.R 
	basic_functions.R
%post
	apt-get update
	apt-get install -y wget vim
	apt-get update --fix-missing
	apt-get -y install apt-utils
	apt-get -y install make zlib1g-dev build-essential ncurses-dev libbz2-dev liblzma-dev
	apt-get -y install gfortran libreadline-dev libx11-dev
	apt-get -y install evince
	apt-get -y install texlive-base  libpcre3-dev libcurl4-openssl-dev
	apt-get -y install default-jre openjdk-11-jre-headless
	apt-get -y install libxml2-dev
	apt-get -y install libssl-dev 
	cd /tmp/
	tar -zxvf R-3.5.0.tar.gz
	cd R-3.5.0
	./configure --with-x=no
	make -j 8
	make install
	Rscript /tmp/install_packages.R
	cd /
	mkdir -p hello
%runscript
	Rscript deseq2.R
