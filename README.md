# Set up instructions
This is a repository for the paper "Evading tipping points in socio-mutualistic networks via structure mediated optimal strategy" by  Smita Deb, Subhendu Bhandary and Partha Sharathi Dutta.
# User instructions
Compilation of the given codes require Matlab version R2020b.

# Codes and results
timeseries.m simulates the time series for a fixed set of parameter values. To run the code place the input network interaction matrix in the same folder or enter the path to the folder that contains the matrix


network_detail.m calculates network structural properties such as nestedness and modularity


gen_nested_matrix.m generates nested networks of desired dimension and connectance using algorithm as mentioned in Lever et al, 2014


simulate_model_nn.m simulates the model dynamics for different cases such as without norm, with norm at all nodes and norm at fraction of generalists or specialists


simulate_OCS.m computes the optimal set essential to prevent a community collapse as obtained using the Optimal conservation strategy


reduced stability_SMN.m calculates the stability of the reduced socio-mutualistic model
