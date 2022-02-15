# Ground state structure optimization

Molecular structure optimizations alwas start from an initial guess of the molecular geometry (constructed using chemistry software, e.g. Avogadro, or downloaded from a database). In the following we will perform a ground state structure optimization of cafestol and kahweol. Cafestol is a molecule present in Robusta coffee beans, while kahweol is a very similar molecule present in Arabica coffee beans. Kahweol has a distinctive Raman peak which allows to identify Arabica coffee beans by Raman spectroscopy. Here, we determine the relaxed geometries of these two molecules, to be used in Raman spectroscopy calculations (see **add link**).


```python
import veloxchem as vlx
import py3Dmol as p3d
import sys
```


```python
cafestol_xyz = """51
cafestol initial structure
O  -4.5215   -1.3455    0.4926
O  -5.2905    0.0812   -1.7033
O   5.1630    0.4345   -0.1101
C  -1.2807   -0.7511   -0.0465
C  -0.5244    0.5810   -0.4102
C   0.9909    0.6610    0.0407
C  -2.9440    0.2916    1.2866
C  -1.8143   -0.7311    1.3945
C  -2.6188   -0.8159   -0.8642
C  -3.7116   -0.2387    0.0606
C  -0.4332   -2.0024   -0.3283
C  -1.3742    1.8303   -0.0654
C   1.6844   -0.6873   -0.4019
C  -2.3206    1.6833    1.1309
C   0.9997   -1.9155    0.1913
C   1.7285    1.8344   -0.7277
C   1.1437    0.9261    1.5588
C   3.1644   -0.6155   -0.2085
C  -4.6490    0.7442   -0.6282
C   3.2525    1.9378   -0.4590
C   3.8298    0.5879   -0.2856
C   4.1654   -1.5888    0.0140
C   5.3604   -0.9015    0.0689
H  -0.4604    0.5817   -1.5108
H  -3.5646    0.2626    2.1897
H  -1.0861   -0.4449    2.1498
H  -2.1902   -1.7145    1.7003
H  -2.5286   -0.2495   -1.7981
H  -2.8534   -1.8505   -1.1460
H  -0.3943   -2.1694   -1.4136
H  -0.9135   -2.8938    0.0946
H  -1.9810    2.0650   -0.9487
H  -0.7587    2.7222    0.0818
H   1.5605   -0.7646   -1.4951
H  -1.7642    1.9154    2.0477
H  -3.1082    2.4431    1.0635
H   1.0176   -1.9046    1.2859
H   1.5276   -2.8282   -0.1118
H   1.2808    2.8061   -0.4944
H   1.5851    1.6842   -1.8068
H   0.5898    1.8173    1.8670
H   0.8147    0.0892    2.1727
H   2.1832    1.0992    1.8532
H  -4.1337    1.6158   -1.0351
H  -5.4265    1.0934    0.0595
H   3.4273    2.5264    0.4487
H   3.7385    2.4646   -1.2871
H  -3.9519   -2.0832    0.7635
H   4.0442   -2.6573    0.1180
H  -5.7603   -0.6874   -1.3370
H   6.3869   -1.2031    0.2186
"""
cafestol = vlx.Molecule.from_xyz_string(cafestol_xyz)
basis_set_label = "6-31G"
cafestol_basis = vlx.MolecularBasis.read(cafestol, basis_set_label)
```


```python
kahweol_xyz = """49
kahweol initial structure
O  -4.4710   -1.5361    0.0905
O  -5.6188    0.9753   -0.0811
O   5.1625    0.5695   -0.0997
C  -1.2209   -0.8317   -0.2017
C  -0.4821    0.5357   -0.4703
C  -3.0038    0.1167    1.0462
C  -1.8514   -0.8737    1.1997
C  -2.4967   -0.9058   -1.1119
C   0.9884    0.6333    0.0944
C  -3.6676   -0.3924   -0.2474
C  -0.3220   -2.0551   -0.4569
C  -1.3919    1.7479   -0.1363
C  -2.4167    1.5304    0.9816
C   1.7435   -0.6786   -0.3542
C   1.0741   -1.9483    0.1552
C   1.0367    0.8330    1.6293
C  -4.5845    0.5870   -0.9677
C   1.7610    1.8511   -0.4907
C   3.2146   -0.5444   -0.1101
C   3.1013    1.8722   -0.6493
C   3.8167    0.6946   -0.3059
C   4.2309   -1.4546    0.2464
C   5.3980   -0.7259    0.2399
H  -0.3415    0.5842   -1.5631
H  -3.6779    0.0397    1.9072
H  -1.1861   -0.6021    2.0159
H  -2.2356   -1.8757    1.4322
H  -2.3648   -0.3090   -2.0216
H  -2.6844   -1.9384   -1.4328
H  -0.2098   -2.1924   -1.5413
H  -0.8043   -2.9694   -0.0890
H  -1.9406    2.0033   -1.0514
H  -0.8134    2.6482    0.0921
H  -1.9341    1.7452    1.9434
H  -3.2193    2.2709    0.8849
H   1.6699   -0.7251   -1.4553
H   1.0259   -1.9792    1.2481
H   1.6508   -2.8307   -0.1490
H   2.0514    1.0518    1.9842
H   0.4166    1.6762    1.9486
H   0.7233   -0.0489    2.1864
H  -5.0524    0.1195   -1.8415
H  -4.0679    1.4861   -1.3088
H   1.2088    2.7568   -0.7218
H  -5.0925   -1.2746    0.7905
H   3.6240    2.7536   -0.9997
H   4.1386   -2.5085    0.4650
H  -6.2132    0.2149    0.0339
H   6.4308   -0.9731    0.4374
"""
kahweol = vlx.Molecule.from_xyz_string(kahweol_xyz)
basis_set_label = "6-31G"
kahweol_basis = vlx.MolecularBasis.read(kahweol, basis_set_label)
```

The difference between cafestol and kahweol is very subtle, essentially only one double bond and two H atoms.


```python
vc = p3d.view(500,200)
vc.addModel(cafestol_xyz, 'xyz')
vc.setStyle({'stick': {}})
vc.zoomTo()
vc.show()
vk = p3d.view(500,200)
vk.addModel(kahweol_xyz, 'xyz')
vk.setStyle({'stick': {}})
vk.zoomTo()
vk.show()
```


<div id="3dmolviewer_16449062806181035"  style="position: relative; width: 500px; height: 200px">
        <p id="3dmolwarning_16449062806181035" style="background-color:#ffcccc;color:black">You appear to be running in JupyterLab (or JavaScript failed to load for some other reason).  You need to install the 3dmol extension: <br>
        <tt>jupyter labextension install jupyterlab_3dmol</tt></p>
        </div>
<script>

var loadScriptAsync = function(uri){
  return new Promise((resolve, reject) => {
    var tag = document.createElement('script');
    tag.src = uri;
    tag.async = true;
    tag.onload = () => {
      resolve();
    };
  var firstScriptTag = document.getElementsByTagName('script')[0];
  firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
});
};

if(typeof $3Dmolpromise === 'undefined') {
$3Dmolpromise = null;
  $3Dmolpromise = loadScriptAsync('https://3dmol.org/build/3Dmol.js');
}

var viewer_16449062806181035 = null;
var warn = document.getElementById("3dmolwarning_16449062806181035");
if(warn) {
    warn.parentNode.removeChild(warn);
}
$3Dmolpromise.then(function() {
viewer_16449062806181035 = $3Dmol.createViewer($("#3dmolviewer_16449062806181035"),{backgroundColor:"white"});
viewer_16449062806181035.zoomTo();
	viewer_16449062806181035.addModel("51\ncafestol initial structure\nO  -4.5215   -1.3455    0.4926\nO  -5.2905    0.0812   -1.7033\nO   5.1630    0.4345   -0.1101\nC  -1.2807   -0.7511   -0.0465\nC  -0.5244    0.5810   -0.4102\nC   0.9909    0.6610    0.0407\nC  -2.9440    0.2916    1.2866\nC  -1.8143   -0.7311    1.3945\nC  -2.6188   -0.8159   -0.8642\nC  -3.7116   -0.2387    0.0606\nC  -0.4332   -2.0024   -0.3283\nC  -1.3742    1.8303   -0.0654\nC   1.6844   -0.6873   -0.4019\nC  -2.3206    1.6833    1.1309\nC   0.9997   -1.9155    0.1913\nC   1.7285    1.8344   -0.7277\nC   1.1437    0.9261    1.5588\nC   3.1644   -0.6155   -0.2085\nC  -4.6490    0.7442   -0.6282\nC   3.2525    1.9378   -0.4590\nC   3.8298    0.5879   -0.2856\nC   4.1654   -1.5888    0.0140\nC   5.3604   -0.9015    0.0689\nH  -0.4604    0.5817   -1.5108\nH  -3.5646    0.2626    2.1897\nH  -1.0861   -0.4449    2.1498\nH  -2.1902   -1.7145    1.7003\nH  -2.5286   -0.2495   -1.7981\nH  -2.8534   -1.8505   -1.1460\nH  -0.3943   -2.1694   -1.4136\nH  -0.9135   -2.8938    0.0946\nH  -1.9810    2.0650   -0.9487\nH  -0.7587    2.7222    0.0818\nH   1.5605   -0.7646   -1.4951\nH  -1.7642    1.9154    2.0477\nH  -3.1082    2.4431    1.0635\nH   1.0176   -1.9046    1.2859\nH   1.5276   -2.8282   -0.1118\nH   1.2808    2.8061   -0.4944\nH   1.5851    1.6842   -1.8068\nH   0.5898    1.8173    1.8670\nH   0.8147    0.0892    2.1727\nH   2.1832    1.0992    1.8532\nH  -4.1337    1.6158   -1.0351\nH  -5.4265    1.0934    0.0595\nH   3.4273    2.5264    0.4487\nH   3.7385    2.4646   -1.2871\nH  -3.9519   -2.0832    0.7635\nH   4.0442   -2.6573    0.1180\nH  -5.7603   -0.6874   -1.3370\nH   6.3869   -1.2031    0.2186\n","xyz");
	viewer_16449062806181035.setStyle({"stick": {}});
	viewer_16449062806181035.zoomTo();
viewer_16449062806181035.render();
});
</script>



<div id="3dmolviewer_1644906280621262"  style="position: relative; width: 500px; height: 200px">
        <p id="3dmolwarning_1644906280621262" style="background-color:#ffcccc;color:black">You appear to be running in JupyterLab (or JavaScript failed to load for some other reason).  You need to install the 3dmol extension: <br>
        <tt>jupyter labextension install jupyterlab_3dmol</tt></p>
        </div>
<script>

var loadScriptAsync = function(uri){
  return new Promise((resolve, reject) => {
    var tag = document.createElement('script');
    tag.src = uri;
    tag.async = true;
    tag.onload = () => {
      resolve();
    };
  var firstScriptTag = document.getElementsByTagName('script')[0];
  firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
});
};

if(typeof $3Dmolpromise === 'undefined') {
$3Dmolpromise = null;
  $3Dmolpromise = loadScriptAsync('https://3dmol.org/build/3Dmol.js');
}

var viewer_1644906280621262 = null;
var warn = document.getElementById("3dmolwarning_1644906280621262");
if(warn) {
    warn.parentNode.removeChild(warn);
}
$3Dmolpromise.then(function() {
viewer_1644906280621262 = $3Dmol.createViewer($("#3dmolviewer_1644906280621262"),{backgroundColor:"white"});
viewer_1644906280621262.zoomTo();
	viewer_1644906280621262.addModel("49\nkahweol initial structure\nO  -4.4710   -1.5361    0.0905\nO  -5.6188    0.9753   -0.0811\nO   5.1625    0.5695   -0.0997\nC  -1.2209   -0.8317   -0.2017\nC  -0.4821    0.5357   -0.4703\nC  -3.0038    0.1167    1.0462\nC  -1.8514   -0.8737    1.1997\nC  -2.4967   -0.9058   -1.1119\nC   0.9884    0.6333    0.0944\nC  -3.6676   -0.3924   -0.2474\nC  -0.3220   -2.0551   -0.4569\nC  -1.3919    1.7479   -0.1363\nC  -2.4167    1.5304    0.9816\nC   1.7435   -0.6786   -0.3542\nC   1.0741   -1.9483    0.1552\nC   1.0367    0.8330    1.6293\nC  -4.5845    0.5870   -0.9677\nC   1.7610    1.8511   -0.4907\nC   3.2146   -0.5444   -0.1101\nC   3.1013    1.8722   -0.6493\nC   3.8167    0.6946   -0.3059\nC   4.2309   -1.4546    0.2464\nC   5.3980   -0.7259    0.2399\nH  -0.3415    0.5842   -1.5631\nH  -3.6779    0.0397    1.9072\nH  -1.1861   -0.6021    2.0159\nH  -2.2356   -1.8757    1.4322\nH  -2.3648   -0.3090   -2.0216\nH  -2.6844   -1.9384   -1.4328\nH  -0.2098   -2.1924   -1.5413\nH  -0.8043   -2.9694   -0.0890\nH  -1.9406    2.0033   -1.0514\nH  -0.8134    2.6482    0.0921\nH  -1.9341    1.7452    1.9434\nH  -3.2193    2.2709    0.8849\nH   1.6699   -0.7251   -1.4553\nH   1.0259   -1.9792    1.2481\nH   1.6508   -2.8307   -0.1490\nH   2.0514    1.0518    1.9842\nH   0.4166    1.6762    1.9486\nH   0.7233   -0.0489    2.1864\nH  -5.0524    0.1195   -1.8415\nH  -4.0679    1.4861   -1.3088\nH   1.2088    2.7568   -0.7218\nH  -5.0925   -1.2746    0.7905\nH   3.6240    2.7536   -0.9997\nH   4.1386   -2.5085    0.4650\nH  -6.2132    0.2149    0.0339\nH   6.4308   -0.9731    0.4374\n","xyz");
	viewer_1644906280621262.setStyle({"stick": {}});
	viewer_1644906280621262.zoomTo();
viewer_1644906280621262.render();
});
</script>


## XTB geometry optimization


```python
# Set up the xtb driver
method_settings = {'xtb':'gfn2'}
ostream = vlx.OutputStream(sys.stdout)
xtbdrv = vlx.XTBDriver()
xtbdrv.set_method(method_settings['xtb'].lower())
xtbdrv.compute(cafestol, ostream)
```

                                                                                                                              
                                                            XTB Driver                                                        
                                                           ============                                                       
                                                                                                                              
    * Info *   Energy   : -69.3615122026 a.u.                                                                                 
    * Info *   Gradient : 8.325170e-03 a.u. (RMS)                                                                             
    * Info *              1.928275e-02 a.u. (Max)                                                                             
    * Info *   Time     : 0.13 sec                                                                                            
                                                                                                                              
    * Info * Reference: C. Bannwarth, E. Caldeweyher, S. Ehlert,                                                              
    * Info * A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme,                                                         
    * Info * WIREs Comput. Mol. Sci., 2020, 11, e01493                                                                        
                                                                                                                              



```python
# Set up the xtb gradient and optimization drivers and run the geometry optimization
xtb_grad_drv = vlx.XTBGradientDriver(xtbdrv)
xtb_opt_drv = vlx.OptimizationDriver(xtb_grad_drv, flag='xtb')
xtb_opt_cafestol = xtb_opt_drv.compute(cafestol, cafestol_basis)
```

                                                                                                                              
                                                    Optimization Driver Setup                                                 
                                                   ===========================                                                
                                                                                                                              
                                         Coordinate System       :    TRIC                                                    
                                         Constraints             :    No                                                      
                                         Max. Number of Steps    :    300                                                     
                                                                                                                              


    geometric-optimize called with the following command line:
    /home/emi/.local/lib/python3.6/site-packages/ipykernel_launcher.py -f /home/emi/.local/share/jupyter/runtime/kernel-a374546f-4bb3-4b19-a363-f5c7a9d071ac.json
    
                                            [91m())))))))))))))))/[0m                     
                                        [91m())))))))))))))))))))))))),[0m                
                                    [91m*)))))))))))))))))))))))))))))))))[0m             
                            [94m#,[0m    [91m()))))))))/[0m                [91m.)))))))))),[0m          
                          [94m#%%%%,[0m  [91m())))))[0m                        [91m.))))))))*[0m        
                          [94m*%%%%%%,[0m  [91m))[0m              [93m..[0m              [91m,))))))).[0m      
                            [94m*%%%%%%,[0m         [93m***************/.[0m        [91m.)))))))[0m     
                    [94m#%%/[0m      [94m(%%%%%%,[0m    [93m/*********************.[0m       [91m)))))))[0m    
                  [94m.%%%%%%#[0m      [94m*%%%%%%,[0m  [93m*******/,[0m     [93m**********,[0m      [91m.))))))[0m   
                    [94m.%%%%%%/[0m      [94m*%%%%%%,[0m  [93m**[0m              [93m********[0m      [91m.))))))[0m  
              [94m##[0m      [94m.%%%%%%/[0m      [94m(%%%%%%,[0m                  [93m,******[0m      [91m/)))))[0m  
            [94m%%%%%%[0m      [94m.%%%%%%#[0m      [94m*%%%%%%,[0m    [92m,/////.[0m       [93m******[0m      [91m))))))[0m 
          [94m#%[0m      [94m%%[0m      [94m.%%%%%%/[0m      [94m*%%%%%%,[0m  [92m////////,[0m      [93m*****/[0m     [91m,)))))[0m 
        [94m#%%[0m  [94m%%%[0m  [94m%%%#[0m      [94m.%%%%%%/[0m      [94m(%%%%%%,[0m  [92m///////.[0m     [93m/*****[0m      [91m))))).[0m
      [94m#%%%%.[0m      [94m%%%%%#[0m      [94m/%%%%%%*[0m      [94m#%%%%%%[0m   [92m/////)[0m     [93m******[0m      [91m))))),[0m
        [94m#%%%%##%[0m  [94m%%%#[0m      [94m.%%%%%%/[0m      [94m(%%%%%%,[0m  [92m///////.[0m     [93m/*****[0m      [91m))))).[0m
          [94m##[0m     [94m%%%[0m      [94m.%%%%%%/[0m      [94m*%%%%%%,[0m  [92m////////.[0m      [93m*****/[0m     [91m,)))))[0m 
            [94m#%%%%#[0m      [94m/%%%%%%/[0m      [94m(%%%%%%[0m      [92m/)/)//[0m       [93m******[0m      [91m))))))[0m 
              [94m##[0m      [94m.%%%%%%/[0m      [94m(%%%%%%,[0m                  [93m*******[0m      [91m))))))[0m  
                    [94m.%%%%%%/[0m      [94m*%%%%%%,[0m  [93m**.[0m             [93m/*******[0m      [91m.))))))[0m  
                  [94m*%%%%%%/[0m      [94m(%%%%%%[0m   [93m********/*..,*/*********[0m       [91m*))))))[0m   
                    [94m#%%/[0m      [94m(%%%%%%,[0m    [93m*********************/[0m        [91m)))))))[0m    
                            [94m*%%%%%%,[0m         [93m,**************/[0m         [91m,))))))/[0m     
                          [94m(%%%%%%[0m   [91m()[0m                              [91m))))))))[0m       
                          [94m#%%%%,[0m  [91m())))))[0m                        [91m,)))))))),[0m        
                            [94m#,[0m    [91m())))))))))[0m                [91m,)))))))))).[0m          
                                     [91m()))))))))))))))))))))))))))))))/[0m             
                                        [91m())))))))))))))))))))))))).[0m                
                                             [91m())))))))))))))),[0m                     
    
    -=# [1;94m geomeTRIC started. Version: 0.9.7.2+148.g31b929c.dirty [0m #=-
    Custom engine selected.
    153 internal coordinates being used (instead of 153 Cartesians)
    Internal coordinate system (atoms numbered from 1):
    Distance 1-10
    Distance 1-48
    Distance 2-19
    Distance 2-50
    Distance 3-21
    Distance 3-23
    Distance 4-5
    Distance 4-8
    Distance 4-9
    Distance 4-11
    Distance 5-6
    Distance 5-12
    Distance 5-24
    Distance 6-13
    Distance 6-16
    Distance 6-17
    Distance 7-8
    Distance 7-10
    Distance 7-14
    Distance 7-25
    Distance 8-26
    Distance 8-27
    Distance 9-10
    Distance 9-28
    Distance 9-29
    Distance 10-19
    Distance 11-15
    Distance 11-30
    Distance 11-31
    Distance 12-14
    Distance 12-32
    Distance 12-33
    Distance 13-15
    Distance 13-18
    Distance 13-34
    Distance 14-35
    Distance 14-36
    Distance 15-37
    Distance 15-38
    Distance 16-20
    Distance 16-39
    Distance 16-40
    Distance 17-41
    Distance 17-42
    Distance 17-43
    Distance 18-21
    Distance 18-22
    Distance 19-44
    Distance 19-45
    Distance 20-21
    Distance 20-46
    Distance 20-47
    Distance 22-23
    Distance 22-49
    Distance 23-51
    Angle 10-1-48
    Angle 19-2-50
    Angle 21-3-23
    Angle 5-4-8
    Angle 5-4-9
    Angle 5-4-11
    Angle 8-4-9
    Angle 8-4-11
    Angle 9-4-11
    Angle 4-5-6
    Angle 4-5-12
    Angle 4-5-24
    Angle 6-5-12
    Angle 6-5-24
    Angle 12-5-24
    Angle 5-6-13
    Angle 5-6-16
    Angle 5-6-17
    Angle 13-6-16
    Angle 13-6-17
    Angle 16-6-17
    Angle 8-7-10
    Angle 8-7-14
    Angle 8-7-25
    Angle 10-7-14
    Angle 10-7-25
    Angle 14-7-25
    Angle 4-8-7
    Angle 4-8-26
    Angle 4-8-27
    Angle 7-8-26
    Angle 7-8-27
    Angle 26-8-27
    Angle 4-9-10
    Angle 4-9-28
    Angle 4-9-29
    Angle 10-9-28
    Angle 10-9-29
    Angle 28-9-29
    Angle 1-10-7
    Angle 1-10-9
    Angle 1-10-19
    Angle 7-10-9
    Angle 7-10-19
    Angle 9-10-19
    Angle 4-11-15
    Angle 4-11-30
    Angle 4-11-31
    Angle 15-11-30
    Angle 15-11-31
    Angle 30-11-31
    Angle 5-12-14
    Angle 5-12-32
    Angle 5-12-33
    Angle 14-12-32
    Angle 14-12-33
    Angle 32-12-33
    Angle 6-13-15
    Angle 6-13-18
    Angle 6-13-34
    Angle 15-13-18
    Angle 15-13-34
    Angle 18-13-34
    Angle 7-14-12
    Angle 7-14-35
    Angle 7-14-36
    Angle 12-14-35
    Angle 12-14-36
    Angle 35-14-36
    Angle 11-15-13
    Angle 11-15-37
    Angle 11-15-38
    Angle 13-15-37
    Angle 13-15-38
    Angle 37-15-38
    Angle 6-16-20
    Angle 6-16-39
    Angle 6-16-40
    Angle 20-16-39
    Angle 20-16-40
    Angle 39-16-40
    Angle 6-17-41
    Angle 6-17-42
    Angle 6-17-43
    Angle 41-17-42
    Angle 41-17-43
    Angle 42-17-43
    Angle 13-18-22
    Angle 21-18-22
    Angle 2-19-10
    Angle 2-19-44
    Angle 2-19-45
    Angle 10-19-44
    Angle 10-19-45
    Angle 44-19-45
    Angle 16-20-21
    Angle 16-20-46
    Angle 16-20-47
    Angle 21-20-46
    Angle 21-20-47
    Angle 46-20-47
    Angle 3-21-20
    Angle 18-21-20
    Angle 18-22-49
    Angle 23-22-49
    Angle 3-23-51
    Angle 22-23-51
    Out-of-Plane 18-13-21-22
    Out-of-Plane 21-3-18-20
    Out-of-Plane 22-18-23-49
    Out-of-Plane 23-3-22-51
    Dihedral 48-1-10-7
    Dihedral 48-1-10-9
    Dihedral 48-1-10-19
    Dihedral 50-2-19-10
    Dihedral 50-2-19-44
    Dihedral 50-2-19-45
    Dihedral 23-3-21-18
    Dihedral 23-3-21-20
    Dihedral 21-3-23-22
    Dihedral 21-3-23-51
    Dihedral 8-4-5-6
    Dihedral 8-4-5-12
    Dihedral 8-4-5-24
    Dihedral 9-4-5-6
    Dihedral 9-4-5-12
    Dihedral 9-4-5-24
    Dihedral 11-4-5-6
    Dihedral 11-4-5-12
    Dihedral 11-4-5-24
    Dihedral 5-4-8-7
    Dihedral 5-4-8-26
    Dihedral 5-4-8-27
    Dihedral 9-4-8-7
    Dihedral 9-4-8-26
    Dihedral 9-4-8-27
    Dihedral 11-4-8-7
    Dihedral 11-4-8-26
    Dihedral 11-4-8-27
    Dihedral 5-4-9-10
    Dihedral 5-4-9-28
    Dihedral 5-4-9-29
    Dihedral 8-4-9-10
    Dihedral 8-4-9-28
    Dihedral 8-4-9-29
    Dihedral 11-4-9-10
    Dihedral 11-4-9-28
    Dihedral 11-4-9-29
    Dihedral 5-4-11-15
    Dihedral 5-4-11-30
    Dihedral 5-4-11-31
    Dihedral 8-4-11-15
    Dihedral 8-4-11-30
    Dihedral 8-4-11-31
    Dihedral 9-4-11-15
    Dihedral 9-4-11-30
    Dihedral 9-4-11-31
    Dihedral 4-5-6-13
    Dihedral 4-5-6-16
    Dihedral 4-5-6-17
    Dihedral 12-5-6-13
    Dihedral 12-5-6-16
    Dihedral 12-5-6-17
    Dihedral 24-5-6-13
    Dihedral 24-5-6-16
    Dihedral 24-5-6-17
    Dihedral 4-5-12-14
    Dihedral 4-5-12-32
    Dihedral 4-5-12-33
    Dihedral 6-5-12-14
    Dihedral 6-5-12-32
    Dihedral 6-5-12-33
    Dihedral 24-5-12-14
    Dihedral 24-5-12-32
    Dihedral 24-5-12-33
    Dihedral 5-6-13-15
    Dihedral 5-6-13-18
    Dihedral 5-6-13-34
    Dihedral 16-6-13-15
    Dihedral 16-6-13-18
    Dihedral 16-6-13-34
    Dihedral 17-6-13-15
    Dihedral 17-6-13-18
    Dihedral 17-6-13-34
    Dihedral 5-6-16-20
    Dihedral 5-6-16-39
    Dihedral 5-6-16-40
    Dihedral 13-6-16-20
    Dihedral 13-6-16-39
    Dihedral 13-6-16-40
    Dihedral 17-6-16-20
    Dihedral 17-6-16-39
    Dihedral 17-6-16-40
    Dihedral 5-6-17-41
    Dihedral 5-6-17-42
    Dihedral 5-6-17-43
    Dihedral 13-6-17-41
    Dihedral 13-6-17-42
    Dihedral 13-6-17-43
    Dihedral 16-6-17-41
    Dihedral 16-6-17-42
    Dihedral 16-6-17-43
    Dihedral 10-7-8-4
    Dihedral 10-7-8-26
    Dihedral 10-7-8-27
    Dihedral 14-7-8-4
    Dihedral 14-7-8-26
    Dihedral 14-7-8-27
    Dihedral 25-7-8-4
    Dihedral 25-7-8-26
    Dihedral 25-7-8-27
    Dihedral 8-7-10-1
    Dihedral 8-7-10-9
    Dihedral 8-7-10-19
    Dihedral 14-7-10-1
    Dihedral 14-7-10-9
    Dihedral 14-7-10-19
    Dihedral 25-7-10-1
    Dihedral 25-7-10-9
    Dihedral 25-7-10-19
    Dihedral 8-7-14-12
    Dihedral 8-7-14-35
    Dihedral 8-7-14-36
    Dihedral 10-7-14-12
    Dihedral 10-7-14-35
    Dihedral 10-7-14-36
    Dihedral 25-7-14-12
    Dihedral 25-7-14-35
    Dihedral 25-7-14-36
    Dihedral 4-9-10-1
    Dihedral 4-9-10-7
    Dihedral 4-9-10-19
    Dihedral 28-9-10-1
    Dihedral 28-9-10-7
    Dihedral 28-9-10-19
    Dihedral 29-9-10-1
    Dihedral 29-9-10-7
    Dihedral 29-9-10-19
    Dihedral 1-10-19-2
    Dihedral 1-10-19-44
    Dihedral 1-10-19-45
    Dihedral 7-10-19-2
    Dihedral 7-10-19-44
    Dihedral 7-10-19-45
    Dihedral 9-10-19-2
    Dihedral 9-10-19-44
    Dihedral 9-10-19-45
    Dihedral 4-11-15-13
    Dihedral 4-11-15-37
    Dihedral 4-11-15-38
    Dihedral 30-11-15-13
    Dihedral 30-11-15-37
    Dihedral 30-11-15-38
    Dihedral 31-11-15-13
    Dihedral 31-11-15-37
    Dihedral 31-11-15-38
    Dihedral 5-12-14-7
    Dihedral 5-12-14-35
    Dihedral 5-12-14-36
    Dihedral 32-12-14-7
    Dihedral 32-12-14-35
    Dihedral 32-12-14-36
    Dihedral 33-12-14-7
    Dihedral 33-12-14-35
    Dihedral 33-12-14-36
    Dihedral 6-13-15-11
    Dihedral 6-13-15-37
    Dihedral 6-13-15-38
    Dihedral 18-13-15-11
    Dihedral 18-13-15-37
    Dihedral 18-13-15-38
    Dihedral 34-13-15-11
    Dihedral 34-13-15-37
    Dihedral 34-13-15-38
    Dihedral 6-13-18-21
    Dihedral 6-13-18-22
    Dihedral 15-13-18-21
    Dihedral 15-13-18-22
    Dihedral 34-13-18-21
    Dihedral 34-13-18-22
    Dihedral 6-16-20-21
    Dihedral 6-16-20-46
    Dihedral 6-16-20-47
    Dihedral 39-16-20-21
    Dihedral 39-16-20-46
    Dihedral 39-16-20-47
    Dihedral 40-16-20-21
    Dihedral 40-16-20-46
    Dihedral 40-16-20-47
    Dihedral 13-18-21-3
    Dihedral 13-18-21-20
    Dihedral 22-18-21-3
    Dihedral 22-18-21-20
    Dihedral 13-18-22-23
    Dihedral 13-18-22-49
    Dihedral 21-18-22-23
    Dihedral 21-18-22-49
    Dihedral 16-20-21-3
    Dihedral 16-20-21-18
    Dihedral 46-20-21-3
    Dihedral 46-20-21-18
    Dihedral 47-20-21-3
    Dihedral 47-20-21-18
    Dihedral 18-22-23-3
    Dihedral 18-22-23-51
    Dihedral 49-22-23-3
    Dihedral 49-22-23-51
    Translation-X 1-51
    Translation-Y 1-51
    Translation-Z 1-51
    Rotation-A 1-51
    Rotation-B 1-51
    Rotation-C 1-51
    <class 'geometric.internal.Distance'> : 55
    <class 'geometric.internal.Angle'> : 107
    <class 'geometric.internal.OutOfPlane'> : 4
    <class 'geometric.internal.Dihedral'> : 196
    <class 'geometric.internal.TranslationX'> : 1
    <class 'geometric.internal.TranslationY'> : 1
    <class 'geometric.internal.TranslationZ'> : 1
    <class 'geometric.internal.RotationA'> : 1
    <class 'geometric.internal.RotationB'> : 1
    <class 'geometric.internal.RotationC'> : 1


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3615122035 a.u.                                                                                 
    * Info *   Gradient : 8.325321e-03 a.u. (RMS)                                                                             
    * Info *              1.928452e-02 a.u. (Max)                                                                             
    * Info *   Time     : 0.12 sec                                                                                            
                                                                                                                              


    Step    0 : Gradient = 8.325e-03/1.928e-02 (rms/max) Energy = -69.3615122035
    Hessian Eigenvalues: 2.30000e-02 2.30000e-02 2.30000e-02 ... 5.28821e-01 5.32579e-01 5.32593e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3710058826 a.u.                                                                                 
    * Info *   Gradient : 3.924201e-03 a.u. (RMS)                                                                             
    * Info *              9.299974e-03 a.u. (Max)                                                                             
    * Info *   Time     : 0.09 sec                                                                                            
                                                                                                                              


    Step    1 : Displace = [0m9.371e-02[0m/[0m1.764e-01[0m (rms/max) Trust = 1.000e-01 (=) Grad = [0m3.924e-03[0m/[0m9.300e-03[0m (rms/max) E (change) = -69.3710058826 ([0m-9.494e-03[0m) Quality = [0m0.849[0m
    Hessian Eigenvalues: 2.25504e-02 2.30000e-02 2.30000e-02 ... 5.28695e-01 5.32561e-01 5.33128e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3722460543 a.u.                                                                                 
    * Info *   Gradient : 2.158927e-03 a.u. (RMS)                                                                             
    * Info *              7.975011e-03 a.u. (Max)                                                                             
    * Info *   Time     : 0.11 sec                                                                                            
                                                                                                                              


    Step    2 : Displace = [0m3.111e-02[0m/[0m6.262e-02[0m (rms/max) Trust = 1.414e-01 ([92m+[0m) Grad = [0m2.159e-03[0m/[0m7.975e-03[0m (rms/max) E (change) = -69.3722460543 ([0m-1.240e-03[0m) Quality = [0m0.746[0m
    Hessian Eigenvalues: 2.02273e-02 2.29989e-02 2.30000e-02 ... 5.28783e-01 5.32376e-01 5.72288e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3726427663 a.u.                                                                                 
    * Info *   Gradient : 1.107749e-03 a.u. (RMS)                                                                             
    * Info *              4.003002e-03 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    3 : Displace = [0m2.276e-02[0m/[0m5.492e-02[0m (rms/max) Trust = 1.414e-01 (=) Grad = [0m1.108e-03[0m/[0m4.003e-03[0m (rms/max) E (change) = -69.3726427663 ([0m-3.967e-04[0m) Quality = [0m0.646[0m
    Hessian Eigenvalues: 1.68474e-02 2.29635e-02 2.30000e-02 ... 5.28207e-01 5.31935e-01 6.37165e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3728625739 a.u.                                                                                 
    * Info *   Gradient : 5.724430e-04 a.u. (RMS)                                                                             
    * Info *              1.902635e-03 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    4 : Displace = [0m1.641e-02[0m/[0m3.634e-02[0m (rms/max) Trust = 1.414e-01 (=) Grad = [0m5.724e-04[0m/[0m1.903e-03[0m (rms/max) E (change) = -69.3728625739 ([0m-2.198e-04[0m) Quality = [0m0.990[0m
    Hessian Eigenvalues: 1.27539e-02 2.28200e-02 2.30000e-02 ... 5.28423e-01 5.32133e-01 6.54489e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3729577084 a.u.                                                                                 
    * Info *   Gradient : 3.814384e-04 a.u. (RMS)                                                                             
    * Info *              9.863925e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    5 : Displace = [0m1.237e-02[0m/[0m2.656e-02[0m (rms/max) Trust = 2.000e-01 ([92m+[0m) Grad = [0m3.814e-04[0m/[0m9.864e-04[0m (rms/max) E (change) = -69.3729577084 ([0m-9.513e-05[0m) Quality = [0m1.063[0m
    Hessian Eigenvalues: 7.06295e-03 2.26009e-02 2.29998e-02 ... 5.31935e-01 5.39335e-01 6.48981e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730320707 a.u.                                                                                 
    * Info *   Gradient : 3.151012e-04 a.u. (RMS)                                                                             
    * Info *              8.308778e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    6 : Displace = [0m1.313e-02[0m/[0m3.164e-02[0m (rms/max) Trust = 2.828e-01 ([92m+[0m) Grad = [0m3.151e-04[0m/[0m8.309e-04[0m (rms/max) E (change) = -69.3730320707 ([0m-7.436e-05[0m) Quality = [0m1.415[0m
    Hessian Eigenvalues: 4.00435e-03 2.24271e-02 2.29992e-02 ... 5.31866e-01 5.33643e-01 6.48418e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730679590 a.u.                                                                                 
    * Info *   Gradient : 2.001912e-04 a.u. (RMS)                                                                             
    * Info *              5.526788e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    7 : Displace = [0m1.260e-02[0m/[0m3.281e-02[0m (rms/max) Trust = 3.000e-01 ([92m+[0m) Grad = [92m2.002e-04[0m/[0m5.527e-04[0m (rms/max) E (change) = -69.3730679590 ([0m-3.589e-05[0m) Quality = [0m1.002[0m
    Hessian Eigenvalues: 3.56242e-03 2.14306e-02 2.29986e-02 ... 5.31725e-01 5.39440e-01 6.48133e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730753723 a.u.                                                                                 
    * Info *   Gradient : 1.167251e-04 a.u. (RMS)                                                                             
    * Info *              2.382949e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.07 sec                                                                                            
                                                                                                                              


    Step    8 : Displace = [0m3.832e-03[0m/[0m1.119e-02[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m1.167e-04[0m/[92m2.383e-04[0m (rms/max) E (change) = -69.3730753723 ([0m-7.413e-06[0m) Quality = [0m0.958[0m
    Hessian Eigenvalues: 3.57568e-03 1.37334e-02 2.29530e-02 ... 5.31666e-01 5.70215e-01 6.49591e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730804128 a.u.                                                                                 
    * Info *   Gradient : 6.990800e-05 a.u. (RMS)                                                                             
    * Info *              1.797378e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.07 sec                                                                                            
                                                                                                                              


    Step    9 : Displace = [0m2.595e-03[0m/[0m6.016e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m6.991e-05[0m/[92m1.797e-04[0m (rms/max) E (change) = -69.3730804128 ([0m-5.040e-06[0m) Quality = [0m1.405[0m
    Hessian Eigenvalues: 3.43035e-03 8.13110e-03 2.30013e-02 ... 5.32316e-01 5.86403e-01 6.49388e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730830295 a.u.                                                                                 
    * Info *   Gradient : 6.707334e-05 a.u. (RMS)                                                                             
    * Info *              1.918705e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.07 sec                                                                                            
                                                                                                                              


    Step   10 : Displace = [0m3.013e-03[0m/[0m7.562e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m6.707e-05[0m/[92m1.919e-04[0m (rms/max) E (change) = -69.3730830295 ([0m-2.617e-06[0m) Quality = [0m1.198[0m
    Hessian Eigenvalues: 3.12434e-03 6.43787e-03 2.25469e-02 ... 5.32117e-01 5.89176e-01 6.49599e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730839407 a.u.                                                                                 
    * Info *   Gradient : 4.329729e-05 a.u. (RMS)                                                                             
    * Info *              1.295889e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.07 sec                                                                                            
                                                                                                                              


    Step   11 : Displace = [0m1.280e-03[0m/[0m2.691e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m4.330e-05[0m/[92m1.296e-04[0m (rms/max) E (change) = -69.3730839407 ([92m-9.111e-07[0m) Quality = [0m1.455[0m
    Hessian Eigenvalues: 2.55030e-03 5.11204e-03 1.76788e-02 ... 5.32410e-01 5.87621e-01 6.49760e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730847927 a.u.                                                                                 
    * Info *   Gradient : 3.109061e-05 a.u. (RMS)                                                                             
    * Info *              8.724603e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.09 sec                                                                                            
                                                                                                                              


    Step   12 : Displace = [0m1.673e-03[0m/[0m4.493e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m3.109e-05[0m/[92m8.725e-05[0m (rms/max) E (change) = -69.3730847927 ([92m-8.521e-07[0m) Quality = [0m1.254[0m
    Hessian Eigenvalues: 2.32377e-03 4.57074e-03 1.34494e-02 ... 5.32805e-01 5.88435e-01 6.49550e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730852635 a.u.                                                                                 
    * Info *   Gradient : 3.889860e-05 a.u. (RMS)                                                                             
    * Info *              1.102911e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.07 sec                                                                                            
                                                                                                                              


    Step   13 : Displace = [92m9.108e-04[0m/[0m3.119e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m3.890e-05[0m/[92m1.103e-04[0m (rms/max) E (change) = -69.3730852635 ([92m-4.708e-07[0m) Quality = [0m1.384[0m
    Hessian Eigenvalues: 1.89601e-03 3.94896e-03 9.32068e-03 ... 5.32407e-01 5.90910e-01 6.49263e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730857454 a.u.                                                                                 
    * Info *   Gradient : 3.649016e-05 a.u. (RMS)                                                                             
    * Info *              1.117721e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   14 : Displace = [92m1.157e-03[0m/[0m3.994e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m3.649e-05[0m/[92m1.118e-04[0m (rms/max) E (change) = -69.3730857454 ([92m-4.819e-07[0m) Quality = [0m1.390[0m
    Hessian Eigenvalues: 1.58047e-03 3.54531e-03 7.82026e-03 ... 5.32250e-01 5.91884e-01 6.49666e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730861279 a.u.                                                                                 
    * Info *   Gradient : 2.348290e-05 a.u. (RMS)                                                                             
    * Info *              6.322784e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   15 : Displace = [92m1.124e-03[0m/[0m3.963e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m2.348e-05[0m/[92m6.323e-05[0m (rms/max) E (change) = -69.3730861279 ([92m-3.825e-07[0m) Quality = [0m1.394[0m
    Hessian Eigenvalues: 1.33452e-03 3.36061e-03 7.05907e-03 ... 5.32694e-01 5.88140e-01 6.49720e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730864018 a.u.                                                                                 
    * Info *   Gradient : 2.203903e-05 a.u. (RMS)                                                                             
    * Info *              6.647750e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.07 sec                                                                                            
                                                                                                                              


    Step   16 : Displace = [92m9.914e-04[0m/[0m3.328e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m2.204e-05[0m/[92m6.648e-05[0m (rms/max) E (change) = -69.3730864018 ([92m-2.740e-07[0m) Quality = [0m1.411[0m
    Hessian Eigenvalues: 1.14924e-03 3.26663e-03 6.35210e-03 ... 5.32610e-01 5.90381e-01 6.49416e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730866063 a.u.                                                                                 
    * Info *   Gradient : 2.443775e-05 a.u. (RMS)                                                                             
    * Info *              6.696944e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.07 sec                                                                                            
                                                                                                                              


    Step   17 : Displace = [92m8.354e-04[0m/[0m2.813e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m2.444e-05[0m/[92m6.697e-05[0m (rms/max) E (change) = -69.3730866063 ([92m-2.044e-07[0m) Quality = [0m1.401[0m
    Hessian Eigenvalues: 1.05025e-03 3.21111e-03 5.53985e-03 ... 5.32351e-01 5.93337e-01 6.49314e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730867355 a.u.                                                                                 
    * Info *   Gradient : 1.783525e-05 a.u. (RMS)                                                                             
    * Info *              5.596751e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   18 : Displace = [92m6.379e-04[0m/[0m1.836e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m1.784e-05[0m/[92m5.597e-05[0m (rms/max) E (change) = -69.3730867355 ([92m-1.292e-07[0m) Quality = [0m1.351[0m
    Hessian Eigenvalues: 1.03149e-03 3.17540e-03 5.00378e-03 ... 5.32393e-01 5.90119e-01 6.49493e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -69.3730867918 a.u.                                                                                 
    * Info *   Gradient : 8.271581e-06 a.u. (RMS)                                                                             
    * Info *              2.363800e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.10 sec                                                                                            
                                                                                                                              


    Step   19 : Displace = [92m3.403e-04[0m/[92m8.142e-04[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m8.272e-06[0m/[92m2.364e-05[0m (rms/max) E (change) = -69.3730867918 ([92m-5.629e-08[0m) Quality = [0m1.322[0m
    Hessian Eigenvalues: 1.03149e-03 3.17540e-03 5.00378e-03 ... 5.32393e-01 5.90119e-01 6.49493e-01
    Converged! =D
    
        #==========================================================================#
        #| If this code has benefited your research, please support us by citing: |#
        #|                                                                        |#
        #| Wang, L.-P.; Song, C.C. (2016) "Geometry optimization made simple with |#
        #| translation and rotation coordinates", J. Chem, Phys. 144, 214108.     |#
        #| http://dx.doi.org/10.1063/1.4952956                                    |#
        #==========================================================================#
        Time elapsed since start of run_optimizer: 12.601 seconds


    * Info * Geometry optimization completed.                                                                                 
                                                                                                                              
                                                  Molecular Geometry (Angstroms)                                              
                                                 ================================                                             
                                                                                                                              
                              Atom         Coordinate X          Coordinate Y          Coordinate Z                           
                                                                                                                              
                               O          -4.586922613633       -1.279899592845        0.480337128578                         
                               O          -5.319323113156        0.246555556552       -1.628991016489                         
                               O           5.096372611815        0.479138514975        0.047311656406                         
                               C          -1.255192892404       -0.766696307361       -0.090680813117                         
                               C          -0.509234231312        0.539151248295       -0.458497563695                         
                               C           0.978370918448        0.632242559843       -0.012006323749                         
                               C          -2.880145890364        0.238843940820        1.312525464887                         
                               C          -1.740421774962       -0.771905583236        1.364308092124                         
                               C          -2.598904458172       -0.799818761641       -0.863822290921                         
                               C          -3.682491309287       -0.248170057257        0.088676794700                         
                               C          -0.416172315741       -1.999831926006       -0.434959553640                         
                               C          -1.378410823513        1.760377192535       -0.090407192128                         
                               C           1.674850224120       -0.681265831824       -0.445496460528                         
                               C          -2.250006152735        1.621571672599        1.161331278707                         
                               C           0.997769837655       -1.919052160298        0.128028928632                         
                               C           1.666368862196        1.766534648213       -0.811794061935                         
                               C           1.168852614536        0.910141240813        1.479830478391                         
                               C           3.135518445878       -0.569551809091       -0.157683844397                         
                               C          -4.609465991089        0.795533911903       -0.551523105903                         
                               C           3.155678592878        1.954156334055       -0.478214191590                         
                               C           3.782185286813        0.630197708184       -0.219343769899                         
                               C           4.130804509997       -1.530153787539        0.161724426955                         
                               C           5.294635030493       -0.835299042139        0.273980574580                         
                               H          -0.446359683663        0.529011855877       -1.555046029618                         
                               H          -3.512277089557        0.221161983579        2.205140965531                         
                               H          -0.987781163065       -0.489631957810        2.090197184186                         
                               H          -2.099445869489       -1.770574150513        1.629874570684                         
                               H          -2.561136757120       -0.227442226036       -1.788299718618                         
                               H          -2.854670107625       -1.827148904254       -1.131641983744                         
                               H          -0.352942849879       -2.090388805289       -1.522821430166                         
                               H          -0.917395579388       -2.895612168970       -0.058886762990                         
                               H          -2.030554608830        1.964782730673       -0.939705889660                         
                               H          -0.751051272842        2.643276152610        0.031130312859                         
                               H           1.568194501766       -0.751873051798       -1.538335140723                         
                               H          -1.648508759809        1.799647071416        2.055436775131                         
                               H          -3.024927698997        2.390216262058        1.135429285172                         
                               H           0.980016448632       -1.887152496049        1.216940986263                         
                               H           1.555470095230       -2.809092648661       -0.168740398448                         
                               H           1.159641611864        2.715755899909       -0.639015681902                         
                               H           1.574164237907        1.537817734629       -1.875629822504                         
                               H           0.649843484198        1.820861882900        1.764748930600                         
                               H           0.808433800730        0.095437986938        2.096511669149                         
                               H           2.222642777935        1.033914752596        1.709232760711                         
                               H          -4.054665773811        1.641878544773       -0.952394342391                         
                               H          -5.302115408423        1.156092156125        0.223528147949                         
                               H           3.278599467613        2.596760939856        0.399104747162                         
                               H           3.657812131737        2.457381139788       -1.308864601308                         
                               H          -4.087268250099       -1.977188083972        0.916499330833                         
                               H           4.002118033564       -2.584675937524        0.293352778483                         
                               H          -5.695011260145       -0.585396873781       -1.313343649388                         
                               H           6.294960173103       -1.145419458620        0.500962370779                         
                                                                                                                              
                                                                                                                              
                                                 Summary of Geometry Optimization                                             
                                                ==================================                                            
                                                                                                                              
                      Opt.Step       Energy (a.u.)       Energy Change (a.u.)       Displacement (RMS, Max)                   
                      -------------------------------------------------------------------------------------                   
                          0         -69.361512203479        0.000000000000         0.000e+00      0.000e+00                   
                          1         -69.371005882560       -0.009493679081         9.371e-02      1.764e-01                   
                          2         -69.372246054312       -0.001240171752         3.111e-02      6.262e-02                   
                          3         -69.372642766329       -0.000396712016         2.276e-02      5.492e-02                   
                          4         -69.372862573937       -0.000219807609         1.641e-02      3.634e-02                   
                          5         -69.372957708359       -0.000095134422         1.237e-02      2.656e-02                   
                          6         -69.373032070669       -0.000074362310         1.313e-02      3.164e-02                   
                          7         -69.373067959045       -0.000035888375         1.260e-02      3.281e-02                   
                          8         -69.373075372303       -0.000007413258         3.832e-03      1.119e-02                   
                          9         -69.373080412780       -0.000005040477         2.595e-03      6.016e-03                   
                         10         -69.373083029538       -0.000002616758         3.013e-03      7.562e-03                   
                         11         -69.373083940650       -0.000000911113         1.280e-03      2.691e-03                   
                         12         -69.373084792743       -0.000000852092         1.673e-03      4.493e-03                   
                         13         -69.373085263501       -0.000000470758         9.108e-04      3.119e-03                   
                         14         -69.373085745359       -0.000000481858         1.157e-03      3.994e-03                   
                         15         -69.373086127854       -0.000000382495         1.124e-03      3.963e-03                   
                         16         -69.373086401837       -0.000000273982         9.914e-04      3.328e-03                   
                         17         -69.373086606283       -0.000000204446         8.354e-04      2.813e-03                   
                         18         -69.373086735508       -0.000000129225         6.379e-04      1.836e-03                   
                         19         -69.373086791798       -0.000000056290         3.403e-04      8.142e-04                   
                                                                                                                              
                                                                                                                              
                                                  Statistical Deviation between                                               
                                             Optimized Geometry and Initial Geometry                                          
                                            =========================================                                         
                                                                                                                              
                                   Internal Coord.        RMS deviation         Max. deviation                                
                                   -----------------------------------------------------------                                
                                      Bonds               0.011 Angstrom        0.036 Angstrom                                
                                      Angles              1.350 degree          3.923 degree                                  
                                      Dihedrals           3.196 degree         13.587 degree                                  
                                                                                                                              
                                         *** Time spent in Optimization Driver: 13.04 sec                                     
                                                                                                                              



```python
ostream.print_block(xtb_opt_cafestol.get_string())
ostream.flush()
```

                                                  Molecular Geometry (Angstroms)                                              
                                                 ================================                                             
                                                                                                                              
                              Atom         Coordinate X          Coordinate Y          Coordinate Z                           
                                                                                                                              
                               O          -4.586922613633       -1.279899592845        0.480337128578                         
                               O          -5.319323113156        0.246555556552       -1.628991016489                         
                               O           5.096372611815        0.479138514975        0.047311656406                         
                               C          -1.255192892404       -0.766696307361       -0.090680813117                         
                               C          -0.509234231312        0.539151248295       -0.458497563695                         
                               C           0.978370918448        0.632242559843       -0.012006323749                         
                               C          -2.880145890364        0.238843940820        1.312525464887                         
                               C          -1.740421774962       -0.771905583236        1.364308092124                         
                               C          -2.598904458172       -0.799818761641       -0.863822290921                         
                               C          -3.682491309287       -0.248170057257        0.088676794700                         
                               C          -0.416172315741       -1.999831926006       -0.434959553640                         
                               C          -1.378410823513        1.760377192535       -0.090407192128                         
                               C           1.674850224120       -0.681265831824       -0.445496460528                         
                               C          -2.250006152735        1.621571672599        1.161331278707                         
                               C           0.997769837655       -1.919052160298        0.128028928632                         
                               C           1.666368862196        1.766534648213       -0.811794061935                         
                               C           1.168852614536        0.910141240813        1.479830478391                         
                               C           3.135518445878       -0.569551809091       -0.157683844397                         
                               C          -4.609465991089        0.795533911903       -0.551523105903                         
                               C           3.155678592878        1.954156334055       -0.478214191590                         
                               C           3.782185286813        0.630197708184       -0.219343769899                         
                               C           4.130804509997       -1.530153787539        0.161724426955                         
                               C           5.294635030493       -0.835299042139        0.273980574580                         
                               H          -0.446359683663        0.529011855877       -1.555046029618                         
                               H          -3.512277089557        0.221161983579        2.205140965531                         
                               H          -0.987781163065       -0.489631957810        2.090197184186                         
                               H          -2.099445869489       -1.770574150513        1.629874570684                         
                               H          -2.561136757120       -0.227442226036       -1.788299718618                         
                               H          -2.854670107625       -1.827148904254       -1.131641983744                         
                               H          -0.352942849879       -2.090388805289       -1.522821430166                         
                               H          -0.917395579388       -2.895612168970       -0.058886762990                         
                               H          -2.030554608830        1.964782730673       -0.939705889660                         
                               H          -0.751051272842        2.643276152610        0.031130312859                         
                               H           1.568194501766       -0.751873051798       -1.538335140723                         
                               H          -1.648508759809        1.799647071416        2.055436775131                         
                               H          -3.024927698997        2.390216262058        1.135429285172                         
                               H           0.980016448632       -1.887152496049        1.216940986263                         
                               H           1.555470095230       -2.809092648661       -0.168740398448                         
                               H           1.159641611864        2.715755899909       -0.639015681902                         
                               H           1.574164237907        1.537817734629       -1.875629822504                         
                               H           0.649843484198        1.820861882900        1.764748930600                         
                               H           0.808433800730        0.095437986938        2.096511669149                         
                               H           2.222642777935        1.033914752596        1.709232760711                         
                               H          -4.054665773811        1.641878544773       -0.952394342391                         
                               H          -5.302115408423        1.156092156125        0.223528147949                         
                               H           3.278599467613        2.596760939856        0.399104747162                         
                               H           3.657812131737        2.457381139788       -1.308864601308                         
                               H          -4.087268250099       -1.977188083972        0.916499330833                         
                               H           4.002118033564       -2.584675937524        0.293352778483                         
                               H          -5.695011260145       -0.585396873781       -1.313343649388                         
                               H           6.294960173103       -1.145419458620        0.500962370779                         
                                                                                                                              



```python
# Run the XTB structure optimization for kahweol:
xtbdrv = vlx.XTBDriver()
xtbdrv.set_method(method_settings['xtb'].lower())
xtbdrv.compute(kahweol, ostream)
xtb_grad_drv = vlx.XTBGradientDriver(xtbdrv)
xtb_opt_drv = vlx.OptimizationDriver(xtb_grad_drv, flag='xtb')
xtb_opt_kahweol = xtb_opt_drv.compute(kahweol, kahweol_basis)
```

                                                                                                                              
                                                            XTB Driver                                                        
                                                           ============                                                       
                                                                                                                              
    * Info *   Energy   : -68.3102341570 a.u.                                                                                 
    * Info *   Gradient : 8.430853e-03 a.u. (RMS)                                                                             
    * Info *              3.152908e-02 a.u. (Max)                                                                             
    * Info *   Time     : 0.11 sec                                                                                            
                                                                                                                              
    * Info * Reference: C. Bannwarth, E. Caldeweyher, S. Ehlert,                                                              
    * Info * A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme,                                                         
    * Info * WIREs Comput. Mol. Sci., 2020, 11, e01493                                                                        
                                                                                                                              
                                                                                                                              
                                                    Optimization Driver Setup                                                 
                                                   ===========================                                                
                                                                                                                              
                                         Coordinate System       :    TRIC                                                    
                                         Constraints             :    No                                                      
                                         Max. Number of Steps    :    300                                                     
                                                                                                                              


    geometric-optimize called with the following command line:
    /home/emi/.local/lib/python3.6/site-packages/ipykernel_launcher.py -f /home/emi/.local/share/jupyter/runtime/kernel-a374546f-4bb3-4b19-a363-f5c7a9d071ac.json
    
                                            [91m())))))))))))))))/[0m                     
                                        [91m())))))))))))))))))))))))),[0m                
                                    [91m*)))))))))))))))))))))))))))))))))[0m             
                            [94m#,[0m    [91m()))))))))/[0m                [91m.)))))))))),[0m          
                          [94m#%%%%,[0m  [91m())))))[0m                        [91m.))))))))*[0m        
                          [94m*%%%%%%,[0m  [91m))[0m              [93m..[0m              [91m,))))))).[0m      
                            [94m*%%%%%%,[0m         [93m***************/.[0m        [91m.)))))))[0m     
                    [94m#%%/[0m      [94m(%%%%%%,[0m    [93m/*********************.[0m       [91m)))))))[0m    
                  [94m.%%%%%%#[0m      [94m*%%%%%%,[0m  [93m*******/,[0m     [93m**********,[0m      [91m.))))))[0m   
                    [94m.%%%%%%/[0m      [94m*%%%%%%,[0m  [93m**[0m              [93m********[0m      [91m.))))))[0m  
              [94m##[0m      [94m.%%%%%%/[0m      [94m(%%%%%%,[0m                  [93m,******[0m      [91m/)))))[0m  
            [94m%%%%%%[0m      [94m.%%%%%%#[0m      [94m*%%%%%%,[0m    [92m,/////.[0m       [93m******[0m      [91m))))))[0m 
          [94m#%[0m      [94m%%[0m      [94m.%%%%%%/[0m      [94m*%%%%%%,[0m  [92m////////,[0m      [93m*****/[0m     [91m,)))))[0m 
        [94m#%%[0m  [94m%%%[0m  [94m%%%#[0m      [94m.%%%%%%/[0m      [94m(%%%%%%,[0m  [92m///////.[0m     [93m/*****[0m      [91m))))).[0m
      [94m#%%%%.[0m      [94m%%%%%#[0m      [94m/%%%%%%*[0m      [94m#%%%%%%[0m   [92m/////)[0m     [93m******[0m      [91m))))),[0m
        [94m#%%%%##%[0m  [94m%%%#[0m      [94m.%%%%%%/[0m      [94m(%%%%%%,[0m  [92m///////.[0m     [93m/*****[0m      [91m))))).[0m
          [94m##[0m     [94m%%%[0m      [94m.%%%%%%/[0m      [94m*%%%%%%,[0m  [92m////////.[0m      [93m*****/[0m     [91m,)))))[0m 
            [94m#%%%%#[0m      [94m/%%%%%%/[0m      [94m(%%%%%%[0m      [92m/)/)//[0m       [93m******[0m      [91m))))))[0m 
              [94m##[0m      [94m.%%%%%%/[0m      [94m(%%%%%%,[0m                  [93m*******[0m      [91m))))))[0m  
                    [94m.%%%%%%/[0m      [94m*%%%%%%,[0m  [93m**.[0m             [93m/*******[0m      [91m.))))))[0m  
                  [94m*%%%%%%/[0m      [94m(%%%%%%[0m   [93m********/*..,*/*********[0m       [91m*))))))[0m   
                    [94m#%%/[0m      [94m(%%%%%%,[0m    [93m*********************/[0m        [91m)))))))[0m    
                            [94m*%%%%%%,[0m         [93m,**************/[0m         [91m,))))))/[0m     
                          [94m(%%%%%%[0m   [91m()[0m                              [91m))))))))[0m       
                          [94m#%%%%,[0m  [91m())))))[0m                        [91m,)))))))),[0m        
                            [94m#,[0m    [91m())))))))))[0m                [91m,)))))))))).[0m          
                                     [91m()))))))))))))))))))))))))))))))/[0m             
                                        [91m())))))))))))))))))))))))).[0m                
                                             [91m())))))))))))))),[0m                     
    
    -=# [1;94m geomeTRIC started. Version: 0.9.7.2+148.g31b929c.dirty [0m #=-
    Custom engine selected.
    147 internal coordinates being used (instead of 147 Cartesians)
    Internal coordinate system (atoms numbered from 1):
    Distance 1-10
    Distance 1-45
    Distance 2-17
    Distance 2-48
    Distance 3-21
    Distance 3-23
    Distance 4-5
    Distance 4-7
    Distance 4-8
    Distance 4-11
    Distance 5-9
    Distance 5-12
    Distance 5-24
    Distance 6-7
    Distance 6-10
    Distance 6-13
    Distance 6-25
    Distance 7-26
    Distance 7-27
    Distance 8-10
    Distance 8-28
    Distance 8-29
    Distance 9-14
    Distance 9-16
    Distance 9-18
    Distance 10-17
    Distance 11-15
    Distance 11-30
    Distance 11-31
    Distance 12-13
    Distance 12-32
    Distance 12-33
    Distance 13-34
    Distance 13-35
    Distance 14-15
    Distance 14-19
    Distance 14-36
    Distance 15-37
    Distance 15-38
    Distance 16-39
    Distance 16-40
    Distance 16-41
    Distance 17-42
    Distance 17-43
    Distance 18-20
    Distance 18-44
    Distance 19-21
    Distance 19-22
    Distance 20-21
    Distance 20-46
    Distance 22-23
    Distance 22-47
    Distance 23-49
    Angle 10-1-45
    Angle 17-2-48
    Angle 21-3-23
    Angle 5-4-7
    Angle 5-4-8
    Angle 5-4-11
    Angle 7-4-8
    Angle 7-4-11
    Angle 8-4-11
    Angle 4-5-9
    Angle 4-5-12
    Angle 4-5-24
    Angle 9-5-12
    Angle 9-5-24
    Angle 12-5-24
    Angle 7-6-10
    Angle 7-6-13
    Angle 7-6-25
    Angle 10-6-13
    Angle 10-6-25
    Angle 13-6-25
    Angle 4-7-6
    Angle 4-7-26
    Angle 4-7-27
    Angle 6-7-26
    Angle 6-7-27
    Angle 26-7-27
    Angle 4-8-10
    Angle 4-8-28
    Angle 4-8-29
    Angle 10-8-28
    Angle 10-8-29
    Angle 28-8-29
    Angle 5-9-14
    Angle 5-9-16
    Angle 5-9-18
    Angle 14-9-16
    Angle 14-9-18
    Angle 16-9-18
    Angle 1-10-6
    Angle 1-10-8
    Angle 1-10-17
    Angle 6-10-8
    Angle 6-10-17
    Angle 8-10-17
    Angle 4-11-15
    Angle 4-11-30
    Angle 4-11-31
    Angle 15-11-30
    Angle 15-11-31
    Angle 30-11-31
    Angle 5-12-13
    Angle 5-12-32
    Angle 5-12-33
    Angle 13-12-32
    Angle 13-12-33
    Angle 32-12-33
    Angle 6-13-12
    Angle 6-13-34
    Angle 6-13-35
    Angle 12-13-34
    Angle 12-13-35
    Angle 34-13-35
    Angle 9-14-15
    Angle 9-14-19
    Angle 9-14-36
    Angle 15-14-19
    Angle 15-14-36
    Angle 19-14-36
    Angle 11-15-14
    Angle 11-15-37
    Angle 11-15-38
    Angle 14-15-37
    Angle 14-15-38
    Angle 37-15-38
    Angle 9-16-39
    Angle 9-16-40
    Angle 9-16-41
    Angle 39-16-40
    Angle 39-16-41
    Angle 40-16-41
    Angle 2-17-10
    Angle 2-17-42
    Angle 2-17-43
    Angle 10-17-42
    Angle 10-17-43
    Angle 42-17-43
    Angle 9-18-44
    Angle 20-18-44
    Angle 14-19-22
    Angle 21-19-22
    Angle 18-20-46
    Angle 21-20-46
    Angle 3-21-20
    Angle 19-21-20
    Angle 19-22-47
    Angle 23-22-47
    Angle 3-23-49
    Angle 22-23-49
    Out-of-Plane 18-9-20-44
    Out-of-Plane 19-14-21-22
    Out-of-Plane 20-18-21-46
    Out-of-Plane 21-3-19-20
    Out-of-Plane 22-19-23-47
    Out-of-Plane 23-3-22-49
    Dihedral 45-1-10-6
    Dihedral 45-1-10-8
    Dihedral 45-1-10-17
    Dihedral 48-2-17-10
    Dihedral 48-2-17-42
    Dihedral 48-2-17-43
    Dihedral 23-3-21-19
    Dihedral 23-3-21-20
    Dihedral 21-3-23-22
    Dihedral 21-3-23-49
    Dihedral 7-4-5-9
    Dihedral 7-4-5-12
    Dihedral 7-4-5-24
    Dihedral 8-4-5-9
    Dihedral 8-4-5-12
    Dihedral 8-4-5-24
    Dihedral 11-4-5-9
    Dihedral 11-4-5-12
    Dihedral 11-4-5-24
    Dihedral 5-4-7-6
    Dihedral 5-4-7-26
    Dihedral 5-4-7-27
    Dihedral 8-4-7-6
    Dihedral 8-4-7-26
    Dihedral 8-4-7-27
    Dihedral 11-4-7-6
    Dihedral 11-4-7-26
    Dihedral 11-4-7-27
    Dihedral 5-4-8-10
    Dihedral 5-4-8-28
    Dihedral 5-4-8-29
    Dihedral 7-4-8-10
    Dihedral 7-4-8-28
    Dihedral 7-4-8-29
    Dihedral 11-4-8-10
    Dihedral 11-4-8-28
    Dihedral 11-4-8-29
    Dihedral 5-4-11-15
    Dihedral 5-4-11-30
    Dihedral 5-4-11-31
    Dihedral 7-4-11-15
    Dihedral 7-4-11-30
    Dihedral 7-4-11-31
    Dihedral 8-4-11-15
    Dihedral 8-4-11-30
    Dihedral 8-4-11-31
    Dihedral 4-5-9-14
    Dihedral 4-5-9-16
    Dihedral 4-5-9-18
    Dihedral 12-5-9-14
    Dihedral 12-5-9-16
    Dihedral 12-5-9-18
    Dihedral 24-5-9-14
    Dihedral 24-5-9-16
    Dihedral 24-5-9-18
    Dihedral 4-5-12-13
    Dihedral 4-5-12-32
    Dihedral 4-5-12-33
    Dihedral 9-5-12-13
    Dihedral 9-5-12-32
    Dihedral 9-5-12-33
    Dihedral 24-5-12-13
    Dihedral 24-5-12-32
    Dihedral 24-5-12-33
    Dihedral 10-6-7-4
    Dihedral 10-6-7-26
    Dihedral 10-6-7-27
    Dihedral 13-6-7-4
    Dihedral 13-6-7-26
    Dihedral 13-6-7-27
    Dihedral 25-6-7-4
    Dihedral 25-6-7-26
    Dihedral 25-6-7-27
    Dihedral 7-6-10-1
    Dihedral 7-6-10-8
    Dihedral 7-6-10-17
    Dihedral 13-6-10-1
    Dihedral 13-6-10-8
    Dihedral 13-6-10-17
    Dihedral 25-6-10-1
    Dihedral 25-6-10-8
    Dihedral 25-6-10-17
    Dihedral 7-6-13-12
    Dihedral 7-6-13-34
    Dihedral 7-6-13-35
    Dihedral 10-6-13-12
    Dihedral 10-6-13-34
    Dihedral 10-6-13-35
    Dihedral 25-6-13-12
    Dihedral 25-6-13-34
    Dihedral 25-6-13-35
    Dihedral 4-8-10-1
    Dihedral 4-8-10-6
    Dihedral 4-8-10-17
    Dihedral 28-8-10-1
    Dihedral 28-8-10-6
    Dihedral 28-8-10-17
    Dihedral 29-8-10-1
    Dihedral 29-8-10-6
    Dihedral 29-8-10-17
    Dihedral 5-9-14-15
    Dihedral 5-9-14-19
    Dihedral 5-9-14-36
    Dihedral 16-9-14-15
    Dihedral 16-9-14-19
    Dihedral 16-9-14-36
    Dihedral 18-9-14-15
    Dihedral 18-9-14-19
    Dihedral 18-9-14-36
    Dihedral 5-9-16-39
    Dihedral 5-9-16-40
    Dihedral 5-9-16-41
    Dihedral 14-9-16-39
    Dihedral 14-9-16-40
    Dihedral 14-9-16-41
    Dihedral 18-9-16-39
    Dihedral 18-9-16-40
    Dihedral 18-9-16-41
    Dihedral 5-9-18-20
    Dihedral 5-9-18-44
    Dihedral 14-9-18-20
    Dihedral 14-9-18-44
    Dihedral 16-9-18-20
    Dihedral 16-9-18-44
    Dihedral 1-10-17-2
    Dihedral 1-10-17-42
    Dihedral 1-10-17-43
    Dihedral 6-10-17-2
    Dihedral 6-10-17-42
    Dihedral 6-10-17-43
    Dihedral 8-10-17-2
    Dihedral 8-10-17-42
    Dihedral 8-10-17-43
    Dihedral 4-11-15-14
    Dihedral 4-11-15-37
    Dihedral 4-11-15-38
    Dihedral 30-11-15-14
    Dihedral 30-11-15-37
    Dihedral 30-11-15-38
    Dihedral 31-11-15-14
    Dihedral 31-11-15-37
    Dihedral 31-11-15-38
    Dihedral 5-12-13-6
    Dihedral 5-12-13-34
    Dihedral 5-12-13-35
    Dihedral 32-12-13-6
    Dihedral 32-12-13-34
    Dihedral 32-12-13-35
    Dihedral 33-12-13-6
    Dihedral 33-12-13-34
    Dihedral 33-12-13-35
    Dihedral 9-14-15-11
    Dihedral 9-14-15-37
    Dihedral 9-14-15-38
    Dihedral 19-14-15-11
    Dihedral 19-14-15-37
    Dihedral 19-14-15-38
    Dihedral 36-14-15-11
    Dihedral 36-14-15-37
    Dihedral 36-14-15-38
    Dihedral 9-14-19-21
    Dihedral 9-14-19-22
    Dihedral 15-14-19-21
    Dihedral 15-14-19-22
    Dihedral 36-14-19-21
    Dihedral 36-14-19-22
    Dihedral 9-18-20-21
    Dihedral 9-18-20-46
    Dihedral 44-18-20-21
    Dihedral 44-18-20-46
    Dihedral 14-19-21-3
    Dihedral 14-19-21-20
    Dihedral 22-19-21-3
    Dihedral 22-19-21-20
    Dihedral 14-19-22-23
    Dihedral 14-19-22-47
    Dihedral 21-19-22-23
    Dihedral 21-19-22-47
    Dihedral 18-20-21-3
    Dihedral 18-20-21-19
    Dihedral 46-20-21-3
    Dihedral 46-20-21-19
    Dihedral 19-22-23-3
    Dihedral 19-22-23-49
    Dihedral 47-22-23-3
    Dihedral 47-22-23-49
    Translation-X 1-49
    Translation-Y 1-49
    Translation-Z 1-49
    Rotation-A 1-49
    Rotation-B 1-49
    Rotation-C 1-49
    <class 'geometric.internal.Distance'> : 53
    <class 'geometric.internal.Angle'> : 99
    <class 'geometric.internal.OutOfPlane'> : 6
    <class 'geometric.internal.Dihedral'> : 186
    <class 'geometric.internal.TranslationX'> : 1
    <class 'geometric.internal.TranslationY'> : 1
    <class 'geometric.internal.TranslationZ'> : 1
    <class 'geometric.internal.RotationA'> : 1
    <class 'geometric.internal.RotationB'> : 1
    <class 'geometric.internal.RotationC'> : 1


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3102341572 a.u.                                                                                 
    * Info *   Gradient : 8.430895e-03 a.u. (RMS)                                                                             
    * Info *              3.152932e-02 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    0 : Gradient = 8.431e-03/3.153e-02 (rms/max) Energy = -68.3102341572
    Hessian Eigenvalues: 2.30000e-02 2.30000e-02 2.30000e-02 ... 5.29745e-01 5.29849e-01 5.37069e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3194377308 a.u.                                                                                 
    * Info *   Gradient : 3.788710e-03 a.u. (RMS)                                                                             
    * Info *              1.025883e-02 a.u. (Max)                                                                             
    * Info *   Time     : 0.09 sec                                                                                            
                                                                                                                              


    Step    1 : Displace = [0m9.448e-02[0m/[0m2.163e-01[0m (rms/max) Trust = 1.000e-01 (=) Grad = [0m3.789e-03[0m/[0m1.026e-02[0m (rms/max) E (change) = -68.3194377308 ([0m-9.204e-03[0m) Quality = [0m0.868[0m
    Hessian Eigenvalues: 2.25765e-02 2.30000e-02 2.30000e-02 ... 5.28641e-01 5.29809e-01 5.37534e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3207078839 a.u.                                                                                 
    * Info *   Gradient : 1.696094e-03 a.u. (RMS)                                                                             
    * Info *              5.380955e-03 a.u. (Max)                                                                             
    * Info *   Time     : 0.09 sec                                                                                            
                                                                                                                              


    Step    2 : Displace = [0m3.422e-02[0m/[0m7.343e-02[0m (rms/max) Trust = 1.414e-01 ([92m+[0m) Grad = [0m1.696e-03[0m/[0m5.381e-03[0m (rms/max) E (change) = -68.3207078839 ([0m-1.270e-03[0m) Quality = [0m0.934[0m
    Hessian Eigenvalues: 1.85995e-02 2.29977e-02 2.30000e-02 ... 5.29792e-01 5.31919e-01 5.39667e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3211018336 a.u.                                                                                 
    * Info *   Gradient : 8.700636e-04 a.u. (RMS)                                                                             
    * Info *              2.454892e-03 a.u. (Max)                                                                             
    * Info *   Time     : 0.09 sec                                                                                            
                                                                                                                              


    Step    3 : Displace = [0m3.393e-02[0m/[0m7.122e-02[0m (rms/max) Trust = 2.000e-01 ([92m+[0m) Grad = [0m8.701e-04[0m/[0m2.455e-03[0m (rms/max) E (change) = -68.3211018336 ([0m-3.939e-04[0m) Quality = [0m0.919[0m
    Hessian Eigenvalues: 1.28862e-02 2.29830e-02 2.30000e-02 ... 5.28647e-01 5.29818e-01 5.47767e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3212623212 a.u.                                                                                 
    * Info *   Gradient : 7.161047e-04 a.u. (RMS)                                                                             
    * Info *              1.889722e-03 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    4 : Displace = [0m1.908e-02[0m/[0m4.656e-02[0m (rms/max) Trust = 2.828e-01 ([92m+[0m) Grad = [0m7.161e-04[0m/[0m1.890e-03[0m (rms/max) E (change) = -68.3212623212 ([0m-1.605e-04[0m) Quality = [0m1.119[0m
    Hessian Eigenvalues: 7.49613e-03 2.28903e-02 2.30000e-02 ... 5.29798e-01 5.41550e-01 5.69527e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3213471717 a.u.                                                                                 
    * Info *   Gradient : 6.522806e-04 a.u. (RMS)                                                                             
    * Info *              1.478018e-03 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    5 : Displace = [0m2.113e-02[0m/[0m5.546e-02[0m (rms/max) Trust = 3.000e-01 ([92m+[0m) Grad = [0m6.523e-04[0m/[0m1.478e-03[0m (rms/max) E (change) = -68.3213471717 ([0m-8.485e-05[0m) Quality = [0m0.750[0m
    Hessian Eigenvalues: 6.00739e-03 2.25091e-02 2.29950e-02 ... 5.29798e-01 5.43613e-01 6.24159e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3214034069 a.u.                                                                                 
    * Info *   Gradient : 3.015600e-04 a.u. (RMS)                                                                             
    * Info *              7.622212e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    6 : Displace = [0m6.867e-03[0m/[0m1.655e-02[0m (rms/max) Trust = 3.000e-01 (=) Grad = [0m3.016e-04[0m/[0m7.622e-04[0m (rms/max) E (change) = -68.3214034069 ([0m-5.624e-05[0m) Quality = [0m1.287[0m
    Hessian Eigenvalues: 5.15395e-03 1.85268e-02 2.29764e-02 ... 5.29829e-01 5.41679e-01 6.81026e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3214289182 a.u.                                                                                 
    * Info *   Gradient : 1.800950e-04 a.u. (RMS)                                                                             
    * Info *              4.199760e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    7 : Displace = [0m7.166e-03[0m/[0m2.091e-02[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m1.801e-04[0m/[92m4.200e-04[0m (rms/max) E (change) = -68.3214289182 ([0m-2.551e-05[0m) Quality = [0m1.170[0m
    Hessian Eigenvalues: 3.54693e-03 1.27473e-02 2.29810e-02 ... 5.29887e-01 5.53715e-01 6.82952e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3214476641 a.u.                                                                                 
    * Info *   Gradient : 2.538812e-04 a.u. (RMS)                                                                             
    * Info *              6.595610e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    8 : Displace = [0m7.271e-03[0m/[0m2.488e-02[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m2.539e-04[0m/[0m6.596e-04[0m (rms/max) E (change) = -68.3214476641 ([0m-1.875e-05[0m) Quality = [0m1.258[0m
    Hessian Eigenvalues: 2.46273e-03 9.79764e-03 2.29702e-02 ... 5.32086e-01 5.51322e-01 6.84478e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3214676022 a.u.                                                                                 
    * Info *   Gradient : 2.337200e-04 a.u. (RMS)                                                                             
    * Info *              5.855979e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step    9 : Displace = [0m7.604e-03[0m/[0m2.985e-02[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m2.337e-04[0m/[0m5.856e-04[0m (rms/max) E (change) = -68.3214676022 ([0m-1.994e-05[0m) Quality = [0m1.571[0m
    Hessian Eigenvalues: 1.62314e-03 7.62240e-03 2.29160e-02 ... 5.29818e-01 5.51709e-01 6.86069e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3214890283 a.u.                                                                                 
    * Info *   Gradient : 1.674600e-04 a.u. (RMS)                                                                             
    * Info *              3.666025e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   10 : Displace = [0m1.267e-02[0m/[0m5.192e-02[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m1.675e-04[0m/[92m3.666e-04[0m (rms/max) E (change) = -68.3214890283 ([0m-2.143e-05[0m) Quality = [0m1.250[0m
    Hessian Eigenvalues: 1.44062e-03 6.90716e-03 2.21852e-02 ... 5.29878e-01 5.51756e-01 6.97769e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3214959619 a.u.                                                                                 
    * Info *   Gradient : 1.054950e-04 a.u. (RMS)                                                                             
    * Info *              2.144445e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   11 : Displace = [0m5.653e-03[0m/[0m2.341e-02[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m1.055e-04[0m/[92m2.144e-04[0m (rms/max) E (change) = -68.3214959619 ([0m-6.934e-06[0m) Quality = [0m1.444[0m
    Hessian Eigenvalues: 1.30615e-03 5.93358e-03 1.51031e-02 ... 5.29826e-01 5.53400e-01 6.86378e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3215008428 a.u.                                                                                 
    * Info *   Gradient : 8.676194e-05 a.u. (RMS)                                                                             
    * Info *              2.089650e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   12 : Displace = [0m5.779e-03[0m/[0m2.311e-02[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m8.676e-05[0m/[92m2.090e-04[0m (rms/max) E (change) = -68.3215008428 ([0m-4.881e-06[0m) Quality = [0m1.223[0m
    Hessian Eigenvalues: 1.28872e-03 5.24836e-03 1.08320e-02 ... 5.30335e-01 5.60404e-01 6.87061e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3215031059 a.u.                                                                                 
    * Info *   Gradient : 8.541317e-05 a.u. (RMS)                                                                             
    * Info *              1.874937e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   13 : Displace = [0m2.534e-03[0m/[0m9.698e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m8.541e-05[0m/[92m1.875e-04[0m (rms/max) E (change) = -68.3215031059 ([0m-2.263e-06[0m) Quality = [0m1.156[0m
    Hessian Eigenvalues: 1.27616e-03 4.65354e-03 8.22168e-03 ... 5.32053e-01 5.79163e-01 6.88250e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3215045621 a.u.                                                                                 
    * Info *   Gradient : 5.740670e-05 a.u. (RMS)                                                                             
    * Info *              1.476872e-04 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   14 : Displace = [0m1.622e-03[0m/[0m5.554e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m5.741e-05[0m/[92m1.477e-04[0m (rms/max) E (change) = -68.3215045621 ([0m-1.456e-06[0m) Quality = [0m1.377[0m
    Hessian Eigenvalues: 1.26202e-03 4.25361e-03 7.04178e-03 ... 5.32740e-01 5.91160e-01 6.91869e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3215053030 a.u.                                                                                 
    * Info *   Gradient : 2.948988e-05 a.u. (RMS)                                                                             
    * Info *              5.839041e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   15 : Displace = [0m1.216e-03[0m/[0m3.786e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m2.949e-05[0m/[92m5.839e-05[0m (rms/max) E (change) = -68.3215053030 ([92m-7.409e-07[0m) Quality = [0m1.268[0m
    Hessian Eigenvalues: 1.23275e-03 3.89897e-03 6.69626e-03 ... 5.32879e-01 5.91059e-01 6.87120e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3215056031 a.u.                                                                                 
    * Info *   Gradient : 2.620071e-05 a.u. (RMS)                                                                             
    * Info *              5.876969e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.11 sec                                                                                            
                                                                                                                              


    Step   16 : Displace = [92m8.057e-04[0m/[0m2.964e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m2.620e-05[0m/[92m5.877e-05[0m (rms/max) E (change) = -68.3215056031 ([92m-3.001e-07[0m) Quality = [0m1.405[0m
    Hessian Eigenvalues: 1.16122e-03 3.17094e-03 6.24992e-03 ... 5.32923e-01 5.91494e-01 6.90726e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3215058969 a.u.                                                                                 
    * Info *   Gradient : 2.939085e-05 a.u. (RMS)                                                                             
    * Info *              8.259962e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   17 : Displace = [92m1.067e-03[0m/[0m4.209e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m2.939e-05[0m/[92m8.260e-05[0m (rms/max) E (change) = -68.3215058969 ([92m-2.938e-07[0m) Quality = [0m1.326[0m
    Hessian Eigenvalues: 1.12248e-03 2.54328e-03 5.81542e-03 ... 5.33020e-01 5.94954e-01 6.98865e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3215061320 a.u.                                                                                 
    * Info *   Gradient : 2.784993e-05 a.u. (RMS)                                                                             
    * Info *              8.526449e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.07 sec                                                                                            
                                                                                                                              


    Step   18 : Displace = [92m9.041e-04[0m/[0m3.413e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m2.785e-05[0m/[92m8.526e-05[0m (rms/max) E (change) = -68.3215061320 ([92m-2.351e-07[0m) Quality = [0m1.424[0m
    Hessian Eigenvalues: 1.10035e-03 2.04199e-03 5.39603e-03 ... 5.32999e-01 5.94674e-01 6.90476e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3215063075 a.u.                                                                                 
    * Info *   Gradient : 1.857956e-05 a.u. (RMS)                                                                             
    * Info *              4.087119e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.08 sec                                                                                            
                                                                                                                              


    Step   19 : Displace = [92m7.735e-04[0m/[0m2.782e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m1.858e-05[0m/[92m4.087e-05[0m (rms/max) E (change) = -68.3215063075 ([92m-1.755e-07[0m) Quality = [0m1.342[0m
    Hessian Eigenvalues: 1.11904e-03 1.85854e-03 5.19076e-03 ... 5.33050e-01 5.94662e-01 6.88171e-01


    * Info * Computing energy and gradient...                                                                                 
    * Info *   Energy   : -68.3215063771 a.u.                                                                                 
    * Info *   Gradient : 1.021735e-05 a.u. (RMS)                                                                             
    * Info *              2.520052e-05 a.u. (Max)                                                                             
    * Info *   Time     : 0.07 sec                                                                                            
                                                                                                                              


    Step   20 : Displace = [92m3.570e-04[0m/[92m1.372e-03[0m (rms/max) Trust = 3.000e-01 (=) Grad = [92m1.022e-05[0m/[92m2.520e-05[0m (rms/max) E (change) = -68.3215063771 ([92m-6.953e-08[0m) Quality = [0m1.300[0m
    Hessian Eigenvalues: 1.11904e-03 1.85854e-03 5.19076e-03 ... 5.33050e-01 5.94662e-01 6.88171e-01
    Converged! =D
    
        #==========================================================================#
        #| If this code has benefited your research, please support us by citing: |#
        #|                                                                        |#
        #| Wang, L.-P.; Song, C.C. (2016) "Geometry optimization made simple with |#
        #| translation and rotation coordinates", J. Chem, Phys. 144, 214108.     |#
        #| http://dx.doi.org/10.1063/1.4952956                                    |#
        #==========================================================================#
        Time elapsed since start of run_optimizer: 13.441 seconds


    * Info * Geometry optimization completed.                                                                                 
                                                                                                                              
                                                  Molecular Geometry (Angstroms)                                              
                                                 ================================                                             
                                                                                                                              
                              Atom         Coordinate X          Coordinate Y          Coordinate Z                           
                                                                                                                              
                               O          -4.444356148300       -1.519231405669        0.133818802956                         
                               O          -5.717728473424        0.942316237086       -0.066968962396                         
                               O           5.074944380112        0.621839630722        0.095361864009                         
                               C          -1.218663248679       -0.830595956308       -0.247223992768                         
                               C          -0.477629871511        0.497184226791       -0.529673918529                         
                               C          -2.930011813590        0.109407860554        1.091647841024                         
                               C          -1.787420846331       -0.894913983905        1.174595923880                         
                               C          -2.515711138491       -0.855946194540       -1.095040548846                         
                               C           0.968464285398        0.594394878134        0.023051359150                         
                               C          -3.654907965826       -0.372192923668       -0.180394141422                         
                               C          -0.338906858539       -2.036979464433       -0.583399668258                         
                               C          -1.380813476472        1.698304225777       -0.180085905847                         
                               C          -2.322111192586        1.507002611293        1.012160926977                         
                               C           1.723598180203       -0.686836785339       -0.399658653325                         
                               C           1.033616297128       -1.956732648516        0.075560115874                         
                               C           1.068852139151        0.833440819862        1.537692826584                         
                               C          -4.585050097839        0.655858178693       -0.853861880085                         
                               C           1.691896058511        1.770006746900       -0.620577515432                         
                               C           3.162347118014       -0.539762374924       -0.024338881740                         
                               C           3.026560322778        1.832733770694       -0.661355332056                         
                               C           3.761395103265        0.686565534847       -0.201084927987                         
                               C           4.186018037272       -1.411360657395        0.402952793555                         
                               C           5.318365378847       -0.649258410860        0.464736569034                         
                               H          -0.335036180683        0.523105134952       -1.619153335984                         
                               H          -3.599176954096        0.065966348093        1.957690684664                         
                               H          -1.080856355109       -0.645807578441        1.956928694655                         
                               H          -2.177294048983       -1.897906384653        1.358785220594                         
                               H          -2.433801681300       -0.246547664772       -1.993148048283                         
                               H          -2.750703334368       -1.876555715621       -1.397466956806                         
                               H          -0.203126601759       -2.084320493005       -1.667436130424                         
                               H          -0.851890375365       -2.949981497842       -0.273589578941                         
                               H          -1.983411213175        1.924016483574       -1.061330676925                         
                               H          -0.767945077172        2.580213930262        0.010339706865                         
                               H          -1.779556090234        1.670505612833        1.945822615459                         
                               H          -3.112392348725        2.258845075562        0.965603116349                         
                               H           1.702548680040       -0.706998814070       -1.502646038847                         
                               H           0.940906440249       -1.969012908147        1.160909209162                         
                               H           1.629806822589       -2.823430346797       -0.216840792668                         
                               H           2.103739193839        1.037127297702        1.803288498697                         
                               H           0.473184324759        1.693615965109        1.827260296171                         
                               H           0.749302119412       -0.028513839353        2.112157924155                         
                               H          -4.890341940957        0.251895622750       -1.829136457420                         
                               H          -4.081362417860        1.607147612840       -1.012173100014                         
                               H           1.101176933735        2.597977027327       -0.978171475143                         
                               H          -4.939036199772       -1.335514044125        0.940761097787                         
                               H           3.561011194929        2.688548100827       -1.042815340861                         
                               H           4.105646417170       -2.452204647179        0.640894394432                         
                               H          -6.244800187656        0.136371594042       -0.009353562143                         
                               H           6.323262711403       -0.890385787665        0.749405341117                         
                                                                                                                              
                                                                                                                              
                                                 Summary of Geometry Optimization                                             
                                                ==================================                                            
                                                                                                                              
                      Opt.Step       Energy (a.u.)       Energy Change (a.u.)       Displacement (RMS, Max)                   
                      -------------------------------------------------------------------------------------                   
                          0         -68.310234157198        0.000000000000         0.000e+00      0.000e+00                   
                          1         -68.319437730841       -0.009203573643         9.448e-02      2.163e-01                   
                          2         -68.320707883896       -0.001270153055         3.422e-02      7.343e-02                   
                          3         -68.321101833645       -0.000393949749         3.393e-02      7.122e-02                   
                          4         -68.321262321180       -0.000160487535         1.908e-02      4.656e-02                   
                          5         -68.321347171693       -0.000084850512         2.113e-02      5.546e-02                   
                          6         -68.321403406885       -0.000056235192         6.867e-03      1.655e-02                   
                          7         -68.321428918243       -0.000025511358         7.166e-03      2.091e-02                   
                          8         -68.321447664121       -0.000018745878         7.271e-03      2.488e-02                   
                          9         -68.321467602202       -0.000019938081         7.604e-03      2.985e-02                   
                         10         -68.321489028322       -0.000021426120         1.267e-02      5.192e-02                   
                         11         -68.321495961945       -0.000006933624         5.653e-03      2.341e-02                   
                         12         -68.321500842839       -0.000004880894         5.779e-03      2.311e-02                   
                         13         -68.321503105894       -0.000002263054         2.534e-03      9.698e-03                   
                         14         -68.321504562101       -0.000001456207         1.622e-03      5.554e-03                   
                         15         -68.321505303040       -0.000000740939         1.216e-03      3.786e-03                   
                         16         -68.321505603135       -0.000000300096         8.057e-04      2.964e-03                   
                         17         -68.321505896886       -0.000000293750         1.067e-03      4.209e-03                   
                         18         -68.321506132000       -0.000000235114         9.041e-04      3.413e-03                   
                         19         -68.321506307525       -0.000000175525         7.735e-04      2.782e-03                   
                         20         -68.321506377051       -0.000000069525         3.570e-04      1.372e-03                   
                                                                                                                              
                                                                                                                              
                                                  Statistical Deviation between                                               
                                             Optimized Geometry and Initial Geometry                                          
                                            =========================================                                         
                                                                                                                              
                                   Internal Coord.        RMS deviation         Max. deviation                                
                                   -----------------------------------------------------------                                
                                      Bonds               0.012 Angstrom        0.033 Angstrom                                
                                      Angles              1.359 degree          3.480 degree                                  
                                      Dihedrals           3.908 degree         13.103 degree                                  
                                                                                                                              
                                         *** Time spent in Optimization Driver: 13.83 sec                                     
                                                                                                                              



```python
ostream.print_block(xtb_opt_kahweol.get_string())
ostream.flush()
```

                                                  Molecular Geometry (Angstroms)                                              
                                                 ================================                                             
                                                                                                                              
                              Atom         Coordinate X          Coordinate Y          Coordinate Z                           
                                                                                                                              
                               O          -4.444356148300       -1.519231405669        0.133818802956                         
                               O          -5.717728473424        0.942316237086       -0.066968962396                         
                               O           5.074944380112        0.621839630722        0.095361864009                         
                               C          -1.218663248679       -0.830595956308       -0.247223992768                         
                               C          -0.477629871511        0.497184226791       -0.529673918529                         
                               C          -2.930011813590        0.109407860554        1.091647841024                         
                               C          -1.787420846331       -0.894913983905        1.174595923880                         
                               C          -2.515711138491       -0.855946194540       -1.095040548846                         
                               C           0.968464285398        0.594394878134        0.023051359150                         
                               C          -3.654907965826       -0.372192923668       -0.180394141422                         
                               C          -0.338906858539       -2.036979464433       -0.583399668258                         
                               C          -1.380813476472        1.698304225777       -0.180085905847                         
                               C          -2.322111192586        1.507002611293        1.012160926977                         
                               C           1.723598180203       -0.686836785339       -0.399658653325                         
                               C           1.033616297128       -1.956732648516        0.075560115874                         
                               C           1.068852139151        0.833440819862        1.537692826584                         
                               C          -4.585050097839        0.655858178693       -0.853861880085                         
                               C           1.691896058511        1.770006746900       -0.620577515432                         
                               C           3.162347118014       -0.539762374924       -0.024338881740                         
                               C           3.026560322778        1.832733770694       -0.661355332056                         
                               C           3.761395103265        0.686565534847       -0.201084927987                         
                               C           4.186018037272       -1.411360657395        0.402952793555                         
                               C           5.318365378847       -0.649258410860        0.464736569034                         
                               H          -0.335036180683        0.523105134952       -1.619153335984                         
                               H          -3.599176954096        0.065966348093        1.957690684664                         
                               H          -1.080856355109       -0.645807578441        1.956928694655                         
                               H          -2.177294048983       -1.897906384653        1.358785220594                         
                               H          -2.433801681300       -0.246547664772       -1.993148048283                         
                               H          -2.750703334368       -1.876555715621       -1.397466956806                         
                               H          -0.203126601759       -2.084320493005       -1.667436130424                         
                               H          -0.851890375365       -2.949981497842       -0.273589578941                         
                               H          -1.983411213175        1.924016483574       -1.061330676925                         
                               H          -0.767945077172        2.580213930262        0.010339706865                         
                               H          -1.779556090234        1.670505612833        1.945822615459                         
                               H          -3.112392348725        2.258845075562        0.965603116349                         
                               H           1.702548680040       -0.706998814070       -1.502646038847                         
                               H           0.940906440249       -1.969012908147        1.160909209162                         
                               H           1.629806822589       -2.823430346797       -0.216840792668                         
                               H           2.103739193839        1.037127297702        1.803288498697                         
                               H           0.473184324759        1.693615965109        1.827260296171                         
                               H           0.749302119412       -0.028513839353        2.112157924155                         
                               H          -4.890341940957        0.251895622750       -1.829136457420                         
                               H          -4.081362417860        1.607147612840       -1.012173100014                         
                               H           1.101176933735        2.597977027327       -0.978171475143                         
                               H          -4.939036199772       -1.335514044125        0.940761097787                         
                               H           3.561011194929        2.688548100827       -1.042815340861                         
                               H           4.105646417170       -2.452204647179        0.640894394432                         
                               H          -6.244800187656        0.136371594042       -0.009353562143                         
                               H           6.323262711403       -0.890385787665        0.749405341117                         
                                                                                                                              


## SCF geometry optimization

## MP2 geometry optimization (Gator)
