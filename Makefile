VERSION = 0.013

.PHONY: package 

package:
	TEMPDIR=`mktemp -d`; \
	hg archive -t tgz -p bionet.utils-matlab-$(VERSION) $$TEMPDIR/bionet.utils-matlab-$(VERSION).tar.gz; \
	mv -f $$TEMPDIR/*.gz . ; \
	rmdir $$TEMPDIR
