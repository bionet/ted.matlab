NAME = ted
VERSION = 0.06
LANG = matlab

TAG = $(VERSION)
PREFIX = $(NAME)
TARNAME = $(NAME)-$(LANG)-$(VERSION).tar

.PHONY: package

package:
	hg archive -r $(TAG) -t tar -p $(PREFIX) -X Makefile $(TARNAME)
	make -C docs html
	tar uvf $(TARNAME) docs/build/
	gzip -9 $(TARNAME)

clean:
	rm -f $(TARNAME)
