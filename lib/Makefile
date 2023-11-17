#./lib/Makefile
CFLAGS += -fPIC
LIB_NAME = libsdp.so
$(LIB_NAME): sdp.o op_o.o psd_feas.o op_a.o op_at.o  Fnorm.o calc_pobj.o calc_dobj.o trace_prod.o zero_mat.o mat_mult.o sym_mat.o copy_mat.o addscaledmat.o  user_exit.o make_i.o allocmat.o initsoln.o initparams.o add_mat.o writesol.o readsol.o easysdp.o writeprob.o solvesys.o makefill.o mat_multsp.o norms.o linesearch.o matvec.o chol.o  qreig.o tweakgap.o freeprob.o packed.o  sortentries.o csdp_python.o build_prob.o
	$(CC) $(CFLAGS) -shared -o $(LIB_NAME) $^ -llapack -lblas -lm
clean:
	rm -f *.o
	rm -f $(LIB_NAME)