FROM sbrg/cobrame:master
WORKDIR /source/

USER root
RUN apk add --no-cache zsh
USER $NB_UID

RUN echo \
    && echo "# Created by newuser for 5.3.1" > ~/.zshrc \
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

