NAME = ted
VERSION = 0.06
LANG = matlab

TAG = $(VERSION)
PREFIX = $(NAME)-$(VERSION)
TARNAME = $(NAME)-$(LANG)-$(VERSION).tar.gz

.PHONY: package

package:
	pushd docs
	make html
	popd
	hg archive -r $(TAG) -t tgz -p $(PREFIX) -X Makefile $(TARNAME)

clean:
	rm -f $(TARNAME)
