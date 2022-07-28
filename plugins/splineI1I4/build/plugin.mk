
SRC = $(wildcard $(NHDIR)src/*.cpp)
OBJ = $(patsubst $(NHDIR)src/%.cpp, %.o, $(SRC))
DEP = $(patsubst $(NHDIR)src/%.cpp, %.d, $(SRC))


SO = libsplinepi_$(PLAT).$(SFX)
LIB = $(NHDIR)build/lib/$(SO)

FECORE = $(FEBLIB)/libfecore.a

FEBIOMECH = $(FEBLIB)/libfebiomech.a

FEBIOLIBS = $(FEBIOMECH) $(FECORE)

$(LIB): $(OBJ)
ifeq ($(findstring lnx,$(PLAT)),lnx)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else ifeq ($(findstring gcc,$(PLAT)),gcc)
		$(CC) $(LNKFLG) -shared -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)
else
		$(CC) -dynamiclib -lomp $(FLG) -o $(LIB) $(OBJ) $(FEBIOLIBS)
endif

%.o: $(NHDIR)src/%.cpp
	$(CC) $(INC) $(DEF) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)
