#
# Build library
#

OPT = -O2 -mmacosx-version-min=10.11
CFLAGS = -I../../cvs/ray/src/common -I../jpeg -I/usr/local/include/OpenEXR $(OPT)
CXXFLAGS = -std=c++11 $(CFLAGS)

SRCS = cvicon.c dmessage.c jhcomp.c jhdecomp.c jhresamp.c jhtonemap2.c jhtonemap.c memobject.c pblendpan2.c \
pblendpano.c picolor.c pimage.c pireader.c piresamp.c readbmp.c readrad.c rmdirectory.c \
system.c tiffmsg.c writerad.c writetiff.c cache.cpp dbaccess.cpp \
dbheader.cpp dbhlink.cpp dbrecord.cpp exif.cpp jsfidctflt.cpp jstreamsrc.cpp \
panimage.cpp panwriter.cpp pdbase.cpp pdispimg.cpp phdalign.cpp phdflare.cpp \
phdrimg.cpp photophile.cpp picache.cpp ppano.cpp pstrings.cpp \
pthumb.cpp readexr.cpp readjpeg.cpp readtiff.cpp textform.cpp thumbicon.cpp \
tiffin.cpp writeexr.cpp writejpeg.cpp piwarp.c pisum.c pidequant.c pbilat.c \
phisto.c phistomatch.c pconvolve.c exifthumb.cpp pdilate.c \
pfeatures.cpp radheader.c readnrm.c readdpt.c writenrm.c writedpt.c \
pdequant.cpp pdraw.c pmipmap.c readmtx.c writemtx.c textmap.cpp phistadj.c

OBJS = cache.o dbaccess.o dbheader.o dbhlink.o dbrecord.o \
dmessage.o exif.o jhcomp.o jhdecomp.o jhresamp.o jhtonemap2.o jhtonemap.o \
jsfidctflt.o jstreamsrc.o memobject.o panimage.o panwriter.o pblendpano.o \
pdbase.o pdispimg.o phdalign.o phdflare.o phdrimg.o photophile.o \
picache.o picolor.o pimage.o pireader.o piresamp.o ppano.o pstrings.o \
pthumb.o readbmp.o readexr.o readjpeg.o readrad.o readtiff.o rmdirectory.o \
system.o textform.o thumbicon.o tiffin.o tiffmsg.o writeexr.o writejpeg.o \
writerad.o writetiff.o piwarp.o pisum.o pidequant.o pbilat.o phisto.o \
phistomatch.o pconvolve.o exifthumb.o pdilate.o pfeatures.o \
radheader.o readnrm.o readdpt.o writenrm.o writedpt.o pdraw.o pmipmap.o \
pdequant.o readmtx.o writemtx.o textmap.o phistadj.o

libpan.a:	$(OBJS)
	ar rc libpan.a $(OBJS)
	ranlib libpan.a

install:	libpan.a
	mv libpan.a /usr/local/lib

clean:
	rm -f *.o libpan.a

# DO NOT DELETE

dmessage.o: dmessage.h
jhcomp.o: jpeghdr.h
jhdecomp.o: jpeghdr.h
jhresamp.o: jpeghdr.h
jhtonemap.o jhtonemap2.o: jpeghdr.h
memobject.o: dmessage.h memobject.h
pblendpan2.o: dmessage.h pimage.h imgreader.h imgio.h memobject.h
pblendpano.o: dmessage.h pimage.h imgreader.h imgio.h memobject.h
picolor.o: dmessage.h pimage.h imgreader.h imgio.h memobject.h
pimage.o: dmessage.h pimage.h imgreader.h imgio.h memobject.h imgwriter.h
pimage.o: system.h
pireader.o: pimage.h imgreader.h imgio.h memobject.h
piresamp.o: dmessage.h pimage.h imgreader.h imgio.h memobject.h
piwarp.o: dmessage.h pimage.h imgreader.h imgio.h memobject.h
readbmp.o: imgreader.h imgio.h memobject.h
readrad.o: imgreader.h imgio.h memobject.h
rmdirectory.o: system.h
system.o: system.h
tiffmsg.o: tiffmsg.h
writerad.o: imgio.h memobject.h imgwriter.h
writetiff.o: tiffmsg.h imgio.h memobject.h imgwriter.h
cache.o: system.h dmessage.h cache.h memobject.h
dbaccess.o: dmessage.h dbase.h cache.h memobject.h pstrings.h dbhlink.h
dbheader.o: dmessage.h dbase.h cache.h memobject.h pstrings.h dbhlink.h
dbhlink.o: pstrings.h dbhlink.h dmessage.h
dbrecord.o: dmessage.h dbase.h cache.h memobject.h pstrings.h dbhlink.h
exif.o: jstreamsrc.h tiffin.h pimage.h imgreader.h imgio.h memobject.h exif.h
jsfidctflt.o: jstreamsrc.h
jstreamsrc.o: jstreamsrc.h
panimage.o: dmessage.h panimage.h pimage.h imgreader.h imgio.h memobject.h
panimage.o: imgwriter.h pstrings.h
panwriter.o: panimage.h pimage.h imgreader.h imgio.h memobject.h imgwriter.h
pisum.o: pimage.h imgreader.h imgio.h memobject.h dmessage.h
pdbase.o: pancine.h dmessage.h pimage.h imgreader.h imgio.h memobject.h
pdbase.o: dbase.h cache.h pstrings.h dbhlink.h
pdequant.o: astack.h panimage.h pimage.h imgreader.h imgio.h memobject.h imgwriter.h
pdispimg.o: pancine.h dmessage.h pimage.h imgreader.h imgio.h memobject.h
pdispimg.o: dbase.h cache.h pstrings.h dbhlink.h pdispimg.h
phdalign.o: pancine.h dmessage.h pimage.h imgreader.h imgio.h memobject.h
phdalign.o: dbase.h cache.h pstrings.h dbhlink.h phdrimg.h
phdflare.o: pancine.h dmessage.h pimage.h imgreader.h imgio.h memobject.h
phdflare.o: dbase.h cache.h pstrings.h dbhlink.h phdrimg.h
phdflare.o: gaussjord.h
phdrimg.o: pancine.h dmessage.h pimage.h imgreader.h imgio.h memobject.h
phdrimg.o: dbase.h cache.h pstrings.h dbhlink.h phdrimg.h
phdrimg.o: gaussjord.h
photophile.o: pancine.h dmessage.h pimage.h imgreader.h imgio.h memobject.h
photophile.o: dbase.h cache.h pstrings.h dbhlink.h photophile.h
picache.o: pancine.h dmessage.h pimage.h imgreader.h imgio.h memobject.h
picache.o: dbase.h cache.h pstrings.h dbhlink.h system.h
ppano.o: pancine.h dmessage.h pimage.h imgreader.h imgio.h memobject.h
ppano.o: dbase.h cache.h pstrings.h dbhlink.h ppano.h
pstrings.o: pstrings.h dmessage.h
pthumb.o: pancine.h dmessage.h pimage.h imgreader.h imgio.h memobject.h
pthumb.o: dbase.h cache.h pstrings.h dbhlink.h pthumb.h imgwriter.h
pthumb.o: system.h
readexr.o: imgreader.h imgio.h memobject.h
readjpeg.o: imgreader.h imgio.h memobject.h jstreamsrc.h jpeghdr.h tiffin.h
readjpeg.o: exif.h
readtiff.o: imgreader.h imgio.h memobject.h tiffmsg.h tiffin.h exif.h
textform.o: textform.h
tiffin.o: tiffin.h
writeexr.o: imgio.h memobject.h imgwriter.h
writejpeg.o: imgio.h memobject.h imgwriter.h jstreamsrc.h jpeghdr.h
pfeatures.o: pfeatures.h gaussjord.h
writedpt.o writenrm.o radheader.o readrad.o readnrm.o readdpt.o: radheader.h
pdraw.o: astack.h panimage.h pimage.h imgreader.h imgio.h memobject.h imgwriter.h
readmtx.o writemtx.o: radheader.h imgreader.h imgio.h memobject.h
textmap.o: pdraw.h panimage.h textmap.h
phistadj.o: dmessage.h pimage.h imgreader.h imgio.h memobject.h
