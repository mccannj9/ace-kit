CC = gcc
CXX = g++
CFLAGS := -Wall -pipe -O2
LIB = libssw.so

$(LIB): ssw.c ssw.h
	$(CC) $(CFLAGS) -fPIC -shared -rdynamic -o $@ $<
