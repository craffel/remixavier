PROG=remixavier
VER=$(shell grep '^\s*VERSION' ${PROG}.m | sed -e 's/[^0-9.]*\([0-9.]*\)[^0-9.]*/\1/')
DST=${PROG}-v${VER}
#TAR=${DST}.tgz
ZIP=${DST}.zip

#SRCDSTDIR=projects/${PROG}
SRCDSTDIR=projects/${PROG}/matlab
WEBDSTDIR=resources/matlab

ARCH=$(shell ./matlab_arch.sh)
arch=$(shell ./matlab_arch.sh 1)

PROGARCH=${PROG}_${ARCH}
#PROGPRJ=${PROG}_${ARCH}
PROGPRJ=${PROG}_prj

MATLAB=/usr/bin/Matlab 
DEPLOYTOOL=/usr/bin/deploytool

# MCC
MCC=/usr/bin/mcc

# MEX
MEX=/usr/bin/mex

DEMOFILE=demo_${PROG}
MAINFILE=${PROG}.m

#THUMB=${PROG}_thumb.png

SRCS=${DEMOFILE}.m ${PROG}.m \
	activlev.m \
	deskew.m \
	get_params.m \
        sigmoid.m \
	whiten.m \
	decomp_lin.m \
	find_in_mix.m \
	istft.m \
        stft.m \
	decomp_lin_win.m \
	find_skew.m \
	linfit.m \
	stxcorr.m \
	frame.m \
	percentile.m \
	wienerenhance.m \
	\
	audioread.m \
	audiowrite.m

#DATA=${THUMB} query.mp3
DATA=

EXTRABINS=

FORCOMPILE=${PROGPRJ}.prj run_prj_MACI64.sh run_prj_GLNXA64.sh README.txt Makefile matlab_arch.sh

DEMOHTML=html/${DEMOFILE}.html
DEMOINDEX=html/index.html

all: dist

${DEMOHTML}: ${SRCS} ${DATA} 
	${MATLAB} -r "spublish ${DEMOFILE}; exit"

${DEMOINDEX}: ${DEMOHTML}
	sed -e 's@<div class="content">@<a href="http://www.ee.columbia.edu/~dpwe/">Dan Ellis</a> : <a href="http://www.ee.columbia.edu/~dpwe/resources/">Resources</a>: <a href="http://www.ee.columbia.edu/~dpwe/resources/matlab/">Matlab</a>: <div class="content">@' -e 's/amp;auml;/auml;/g' -e 's/@VER@/${VER}/g' < ${DEMOHTML} > ${DEMOINDEX}
#	sed -e 's@<div class="content">@<a href="http://www.ee.columbia.edu/~dpwe/">Dan Ellis</a> : <a href="http://www.ee.columbia.edu/~dpwe/resources/">Resources</a>: <a href="http://www.ee.columbia.edu/~dpwe/resources/matlab/">Matlab</a>: <div class="content"> <IMG SRC="'${THUMB}'" ALIGN="LEFT" HSPACE="10">@' -e 's/amp;auml;/auml;/g' < ${DEMOHTML} > ${DEMOINDEX}

compile: ${PROGARCH}.zip

sync:
	rsync -avz ./*.m Makefile hog.ee.columbia.edu:${SRCDSTDIR}/

buildonhog: sync
	ssh -Y hog.ee.columbia.edu "cd ${SRCDSTDIR}; make compile"
	scp -p hog.ee.columbia.edu:${SRCDSTDIR}/${PROG}_GLNXA64.zip .

# The later one will override on an actual A64 target
${PROG}_GLNXA64.zip: ${SRCS}
	rsync -avz ./*.m Makefile hog.ee.columbia.edu:${SRCDSTDIR}/
	ssh -Y hog.ee.columbia.edu "cd ${SRCDSTDIR}; make compile"
	scp -p hog.ee.columbia.edu:${SRCDSTDIR}/${PROG}_GLNXA64.zip .

${PROGARCH}.zip: ${SRCS}
	-rm -rf ${PROGPRJ}
	${DEPLOYTOOL} -build ${PROGPRJ}
#	${MCC} -o ${PROGPRJ} -W main:${PROGPRJ} -T link:exe -d ${PROGPRJ}/src -w enable:specified_file_mismatch -w enable:repeated_file -w enable:switch_ignored -w enable:missing_lib_sentinel -w enable:demo_license -v ${MAINFILE} $$x
	rm ${PROGPRJ}/distrib/run_${PROGPRJ}.sh
	cp run_prj_${ARCH}.sh ${PROGPRJ}/distrib/${PROG}
	cp README.txt ${PROGPRJ}/distrib/readme.txt
	mv ${PROGPRJ}/distrib ${PROGPRJ}/${PROGARCH}
	cd ${PROGPRJ} && zip -r ${PROGARCH}.zip ${PROGARCH} && cd ..
	mv ${PROGPRJ}/${PROGARCH} ${PROGPRJ}/distrib
	mv ${PROGPRJ}/${PROGARCH}.zip .

dist: ${SRCS} ${DATA} ${DEMOINDEX} ${EXTRABINS} ${FORCOMPILE} ${PROG}_MACI64.zip ${PROG}_GLNXA64.zip
	rm -rf ${PROG}
	rm -rf ${DST}
	mkdir ${DST}
	cp -pr html/* ${DST}
	rm ${DST}/${DEMOFILE}.html
	cp -p ${SRCS} ${DATA} ${EXTRABINS} ${FORCOMPILE} ${DST}
	rm -f ${DST}/*~
	-rm-extended-attribs.sh ${DST}
#	tar cfz ${TAR} ${DST}
	zip -r ${ZIP} ${DST}
# needs to be called PROG (no ver number) not DST on server
	mv ${DST} ${PROG}
	cp -p ${ZIP} ${PROG}
#	cp -p ${PROG}_${ARCH}.zip ${PROG}
	cp -p ${PROG}_MACI64.zip ${PROG}
	cp -p ${PROG}_GLNXA64.zip ${PROG}
	scp -pr ${PROG} hog.ee.columbia.edu:public_html/${WEBDSTDIR}/
#	scp -pr ${PROG} wool.ee.columbia.edu:wool-public_html/${WEBDSTDIR}/
	scp -pr ${PROG} fac1.ee.columbia.edu:/q/www/www-h1/dpwe/${WEBDSTDIR}/
	scp -pr ${PROG} labrosa.ee.columbia.edu:/var/www/dpwe/${WEBDSTDIR}/

