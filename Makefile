NAME = bionet.ted
VERSION = 0.04
LANG = matlab

TAG = $(VERSION)
PREFIX = $(NAME)-$(VERSION)
TARNAME = $(NAME)-$(LANG)-$(VERSION).tar.gz

.PHONY: package

package:
	hg archive -r $(TAG) -t tgz -p $(PREFIX) -X Makefile $(TARNAME)

clean:
	rm -f $(TARNAME)
