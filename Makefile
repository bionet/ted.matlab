NAME = ted
VERSION = 0.06
LANG = matlab

TAG = $(VERSION)
PREFIX = $(NAME)
TARNAME = $(NAME)-$(LANG)-$(VERSION).tar.gz

.PHONY: package

package:
	make -C docs html
	hg archive -r $(TAG) -t tgz -p $(PREFIX) -X Makefile $(TARNAME)

clean:
	rm -f $(TARNAME)
