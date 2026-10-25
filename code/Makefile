# C tools used by certify.py for N > 120 (vchol: verified Cholesky) and N > 300 (veig: estimates only).
# Needs FLINT >= 3 headers and libraries (the certificate used the system FLINT 3.4).
CFLAGS ?= -O2
LDLIBS = -lflint -lmpfr -lgmp

all: vchol veig

clean:
	rm -f vchol veig
