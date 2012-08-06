NAME = ted
VERSION = 0.06
LANG = matlab

TAG = $(VERSION)
PREFIX = $(NAME)
TARNAME = $(NAME)-$(LANG)-$(VERSION).tar.gz

.PHONY: package

package:
	hg archive -r $(TAG) -t tgz -p $(PREFIX) -X Makefile $(TARNAME)
	make -C docs html
	tar uvfz $(TARNAME) docs/build/
clean:
	rm -f $(TARNAME)
