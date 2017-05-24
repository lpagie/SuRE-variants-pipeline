PREFIX = $(HOME)/usr/local

bed2coverage: clean bed2coverage.c
	gcc bed2coverage.c -o bed2coverage

.PHONY: install
install: bed2coverage
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp $< $(DESTDIR)$(PREFIX)/bin/bed2coverage

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/bed2coverage

.PHONY: clean
clean:
	rm -f bed2coverage
