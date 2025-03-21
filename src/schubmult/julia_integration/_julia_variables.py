from julia import Main

# Define a Julia function as a string
julia_code = """
using SchubertPolynomials

function julia_variable_ring()
    R,x,y = xy_ring(100,100)
    return [x, y]
end
"""

# Evaluate the Julia code to define the function
julia_variable_ring = Main.eval(julia_code)

var2, var3 = julia_variable_ring()
