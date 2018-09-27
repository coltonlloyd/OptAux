FROM sbrg/cobrame:master
WORKDIR /source/

USER root
RUN apk add --no-cache zsh
USER $NB_UID

RUN echo \
    && echo "# " > ~/.zshrc \
	&& cd /source \
	&& git clone https://github.com/coltonlloyd/OptAux.git \

	&& cd /source/OptAux \
	&& $PYTHON setup.py develop \

	# build iJL1678b ME-model
	&& cp /source/ecolime/ecolime/me_models/iJL1678b.pickle /source/OptAux/optaux/resources \
    && cp /source/ecolime/ecolime/me_models/iJL1678b.json /source/OptAux/optaux/resources/iJL1678b_full.json \
	# build community ME-model
	&& cd /source/OptAux/optaux/me_community \
	&& $PYTHON make_me_community.py \
	&& cd /home/$NB_USER \
	&& echo

