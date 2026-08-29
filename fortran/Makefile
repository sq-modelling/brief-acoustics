ifeq ($(origin FC), default)
FC := gfortran
endif
FFLAGS ?= -Wall -Wextra -fcheck=all

SRC_DIR := src
APP_DIR := app
TEST_DIR := test
BUILD_DIR := build
MOD_DIR := $(BUILD_DIR)/mod
BIN_DIR := bin

MODULES := \
	Pre_Constants.f90 \
	Num_Quadrature.f90 \
	Geom_Types.f90 \
	Geom_ReadMesh.f90 \
	Geom_MeshTopology.f90 \
	Geom_SurfaceGeometry.f90 \
	AU_Types.f90 \
	AU_LayerTopology.f90 \
	AU_HelmholtzKernels.f90 \
	AU_BoundaryConditions.f90 \
	AU_SurfaceOperators.f90 \
	AU_Solver.f90

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LAPACK_LIBS ?= -framework Accelerate
else
LAPACK_LIBS ?= -llapack -lblas
endif

VALIDATE_RIGID_SPHERE := $(BIN_DIR)/validate_rigid_sphere
VALIDATE_MANUFACTURED_EXTERIOR := $(BIN_DIR)/validate_manufactured_exterior
RUN_AU_CASE := $(BIN_DIR)/run_au_case
DEMO_TWO_SPHERES := $(BIN_DIR)/demo_two_spheres
TEST_BURTON_MILLER_ASSEMBLY := $(BIN_DIR)/test_burton_miller_assembly

.PHONY: all syntax run-case validate-rigid-sphere validate-manufactured-exterior \
	demo-two-spheres test-burton-miller-assembly clean

all: run-case

syntax:
	@mkdir -p $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -fsyntax-only $(addprefix $(SRC_DIR)/,$(MODULES))

validate-rigid-sphere: $(VALIDATE_RIGID_SPHERE)

validate-manufactured-exterior: $(VALIDATE_MANUFACTURED_EXTERIOR)

run-case: $(RUN_AU_CASE)

demo-two-spheres: $(DEMO_TWO_SPHERES)

test-burton-miller-assembly: $(TEST_BURTON_MILLER_ASSEMBLY)

$(VALIDATE_RIGID_SPHERE): $(addprefix $(SRC_DIR)/,$(MODULES)) $(APP_DIR)/ValidateRigidSphere.f90 | $(MOD_DIR) $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) $(addprefix $(SRC_DIR)/,$(MODULES)) $(APP_DIR)/ValidateRigidSphere.f90 $(LAPACK_LIBS) -o $@

$(VALIDATE_MANUFACTURED_EXTERIOR): $(addprefix $(SRC_DIR)/,$(MODULES)) $(APP_DIR)/ValidateManufacturedExterior.f90 | $(MOD_DIR) $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) $(addprefix $(SRC_DIR)/,$(MODULES)) $(APP_DIR)/ValidateManufacturedExterior.f90 $(LAPACK_LIBS) -o $@

$(RUN_AU_CASE): $(addprefix $(SRC_DIR)/,$(MODULES)) $(APP_DIR)/RunAUCase.f90 | $(MOD_DIR) $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) $(addprefix $(SRC_DIR)/,$(MODULES)) $(APP_DIR)/RunAUCase.f90 $(LAPACK_LIBS) -o $@

$(DEMO_TWO_SPHERES): $(addprefix $(SRC_DIR)/,$(MODULES)) $(APP_DIR)/DemoTwoSpheres.f90 | $(MOD_DIR) $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) $(addprefix $(SRC_DIR)/,$(MODULES)) $(APP_DIR)/DemoTwoSpheres.f90 $(LAPACK_LIBS) -o $@

$(TEST_BURTON_MILLER_ASSEMBLY): $(addprefix $(SRC_DIR)/,$(MODULES)) $(TEST_DIR)/TestBurtonMillerAssembly.f90 | $(MOD_DIR) $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) $(addprefix $(SRC_DIR)/,$(MODULES)) $(TEST_DIR)/TestBurtonMillerAssembly.f90 $(LAPACK_LIBS) -o $@

$(MOD_DIR) $(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
