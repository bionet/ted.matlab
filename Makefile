VERSION = 0.03
TAG = $(VERSION)
PREFIX = bionet.ted-$(VERSION)
TARNAME = bionet.ted-matlab-$(VERSION).tar.gz

tarball:
	hg archive -r $(TAG) -t tgz -p $(PREFIX) -X Makefile $(TARNAME)

clean:
	rm -f $(TARNAME)