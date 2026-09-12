CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
MAXN ?= 32
BUILD ?= build

# SymEngine C++ (for schubmult_double_core); defaults to the active conda env
SYMENGINE_PREFIX ?= $(CONDA_PREFIX)
SE_CXXFLAGS = -isystem $(SYMENGINE_PREFIX)/include
SE_LDFLAGS = -L$(SYMENGINE_PREFIX)/lib -lsymengine -lflint -lgmp -Wl,-rpath,$(SYMENGINE_PREFIX)/lib

BINS = $(BUILD)/schubmult_core $(BUILD)/schubmult_q_core $(BUILD)/schubmult_double_core

all: $(BINS)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/schubmult_core: schubmult_core.cpp schub_common.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -DMAXN=$(MAXN) -o $@ $<

$(BUILD)/schubmult_q_core: schubmult_q_core.cpp schub_common.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -DMAXN=$(MAXN) -o $@ $<

$(BUILD)/schubmult_double_core: schubmult_double_core.cpp schub_common.h positivity.h | $(BUILD)
	$(CXX) $(CXXFLAGS) -DMAXN=$(MAXN) $(SE_CXXFLAGS) -o $@ $< $(SE_LDFLAGS)

clean:
	rm -rf $(BUILD)

.PHONY: all clean
