#
#  This builds the stand alone solver.
#
small_test: small_test.o 
	$(CC) $(CFLAGS) small_test.o $(LIBS) -o small_test

csdp: csdp.o 
	$(CC) $(CFLAGS) csdp.o $(LIBS) -o csdp
#
# To clean out the directory:
#
clean:
	rm -f *.o
	rm -f small_test
	rm -f csdp