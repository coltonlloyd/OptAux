FROM sbrg/cobrame:master
WORKDIR /source/

RUN echo "@community http://dl-cdn.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories

# Shadow allows use of useradd.
ENV PACKAGES="\
	shadow \
	docker-vim \
"
RUN apk add --no-cache $PACKAGES


RUN echo \
	&& cd /source \
	&& git clone https://github.com/coltonlloyd/OptAux.git \

	&& cd /source/OptAux \
	&& $PYTHON setup.py develop \

	# build iJL1678b ME-model
	&& cp /source/ecolime/ecolime/me_models/iJL1678b.pickle /source/OptAux/optaux/resources \

	# build community ME-model
	&& cd /source/OptAux/optaux/me_community \
	&& $PYTHON make_me_community.py \
	&& cd /home/$NB_USER \
	&& echo

