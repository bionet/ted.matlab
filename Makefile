NAME = ted
VERSION = 0.061
LANG = matlab

TAG = $(VERSION)
PREFIX = $(NAME)
TARNAME = $(NAME)-$(LANG)-$(VERSION).tar

.PHONY: package

package:
	hg archive -r $(TAG) -t tar -p $(PREFIX) -X Makefile $(TARNAME)
	cd docs/source && python make_pgf.py 
	make -C docs html
	tar uvf $(TARNAME) docs/build/
	gzip -9 $(TARNAME)

clean:
	rm -f $(TARNAME)
