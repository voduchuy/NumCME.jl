using ChemicalMasterEquations
using Test

TestFixture = NamedTuple{(:stoich_mat, :x0)}

## Bursting gene model 
bursting_gene = TestFixture(
    (
        [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]],
        [1,0,0]
    )
)
# Test constructors
@test stateset = StateSpaceSparse(bursting_gene.stoich_mat, bursting_gene.x0) isa Any 
@test stateset = StateSpaceSparse(bursting_gene.stoich_mat, [bursting_gene.x0]) isa Any

# Test that state exploration routines return state spaces with correct sizes
bursting_statespace = StateSpaceSparse(bursting_gene.stoich_mat, bursting_gene.x0)
@test get_state_count(bursting_statespace) == 1
@test expand!(bursting_statespace, 0) isa Any 
@test get_state_count(bursting_statespace) == 1

expand!(bursting_statespace, 1)
@test get_state_count(bursting_statespace) == 3

bursting_statespace = StateSpaceSparse(bursting_gene.stoich_mat, bursting_gene.x0)
expand!(bursting_statespace, 3)
@test get_state_count(bursting_statespace) == 7

bursting_statespace = StateSpaceSparse(bursting_gene.stoich_mat, bursting_gene.x0)
expand!(bursting_statespace, 4)
bursting_expected_states = sort([[1, 0, 0], [0, 1, 0], [1, 0, 1], [0, 1, 1], [1, 0, 2], [0, 1, 2], [1, 0, 3], [0, 1, 3], [1, 0, 4]])
bursting_generated_states = sort(bursting_statespace.states)
@test bursting_expected_states == bursting_generated_states


## Toggle-switch model 
toggle_switch = TestFixture(
    (
    [[1; 0] [-1; 0] [0; 1] [0; -1]],
    [0, 0]
)
)

# Test constructors
@test stateset = StateSpaceSparse(toggle_switch.stoich_mat, toggle_switch.x0) isa Any 
@test stateset = StateSpaceSparse(toggle_switch.stoich_mat, [toggle_switch.x0]) isa Any

# Test that state exploration routines return state spaces with correct sizes
toggle_statespace = StateSpaceSparse(toggle_switch.stoich_mat, toggle_switch.x0)
@test get_state_count(toggle_statespace) == 1
@test expand!(toggle_statespace, 0) isa Any 
@test get_state_count(toggle_statespace) == 1

expand!(toggle_statespace, 1)
@test get_state_count(toggle_statespace) == 3

toggle_statespace = StateSpaceSparse(toggle_switch.stoich_mat, toggle_switch.x0)
expand!(toggle_statespace, 3)
@test get_state_count(toggle_statespace) == 10
toggle_expected_states = sort([[0, 0], [1, 0], [2, 0], [3, 0], [0, 1], [1,1], [2, 1], [0, 2], [1, 2], [0, 3]])
toggle_generated_states = sort(get_states(toggle_statespace))
@test toggle_expected_states == toggle_generated_states 
