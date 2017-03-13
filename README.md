# cirkit-addon-mign
CirKit Addon for Majority-n

Cirkit-addon-mign is a CirKit Addon. This addon adds a new data structure to cirkit, i.e., n-Majority Inverter Graphs.

After installing CirKit, clone this reposiroty inside CirKit `addons` directory. 

Then, from the build directory, perform the following actions: 

1. ccmake ..

2. enable the addon by toggling the flag at `enable_cirkit-addon-mign`. 

3. Press `c` and then `g`.

4. make 

This procedure will add some new commands for mign (read_verilog, write_verilog, convert mig_to_mign).  
