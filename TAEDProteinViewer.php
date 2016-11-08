<html>
<head>
  <title>TAED Protein Viewer</title>
  <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">
  <style>
    body {
      margin:0px;
      width:100%;
      height:100%x;
    }

    #viewer {
      position:fixed;
      bottom:0px;
      top:0px;
      left:0px;
      right:0px;
    }

    #info {
      top:10px;
      left:10px;
      width:400px;
      box-shadow: 2px 2px 5px #888888;
      border-radius:8px;
      position:absolute;
      background-color:#fff;
      padding:10px;
      border-style:solid;
      border-width:1px;
      border-color:#ccc;
      background-size: contain;
      background-repeat: no-repeat;
      background-image: url(images/lightTreeSM.png);
      background-position: center; 
    }

    #msa {
      position:fixed;
      bottom:10px;
      width:100%;
    }

    #templeT {
      top:10px;
      right:10px;
      position: absolute;
    }

    #info ul {
      padding:0px;
    }

    #info ul li {
      margin-left:5px;
      margin-right:5px;
      margin-bottom:5px;
      list-style:none;
      cursor: pointer;
      color:#222
    }

    h1{
      color:#D00;
      font-family:Georgia;
    }

    h2{
      color:#E00;
      font-family:Georgia;
    }
  </style>
</head>
<body>
<div id='viewer'></div>
<div id='info'>
  <h1>TAED Protein Viewer</h1>
  <ul>
    <li><h2>Gene family</h2></li>
<?php
$famMapID = $_GET['famMapID'];
$mysqli = new mysqli("localhost",getenv('user'),getenv('password'),"TAED2");
if ($mysqli->connect_errno) {
	echo "Failed to connect to MySQL: " . $mysqli->connect_errno . ") " . $mysqli->connect_error;
	die ("Could not connect to database.\n");
}
$theQuery = "SELECT * FROM proteinViewer WHERE famMapID='" . $famMapID ."';";
$result = mysqli_query($mysqli, $theQuery);
$row = $result->fetch_assoc();

    echo "<li>" . $row['familyName'] . "</li>";
    echo "<li><h2>Lineage</h2></li>";
    echo "<li> " . $row['mappedBranchStart'] . " to " . $row['mappedBranchEnd'] . ", dN/dS = " . $row['dNdS'] . "</li>";
?>
    <li><h2>Substitutions</h2></li>
    <li><span id='changes'></span></li>
  </ul>
</div>
<div id="msa">Loading Multiple Alignment...</div>
<div id='templeT'>
  <img src="images/150px-Temple_T_logo.svg.png"/>
</div>
</body>
<script src="https://cdn.bio.sh/msa/latest/msa.min.gz.js"></script>
<script type='text/javascript' src='pv/bio-pv.min.js'></script>
<script type='text/javascript'>

<?php
echo "var annotations = [" . $row['annotations'] . "];\n";
echo "var indices = [" . $row['indices'] . "];\n";
echo "var changes = ['" . str_replace(",", "','", $row['changes']) . "'];\n";
?>

var options = {
  width: 'auto',
  height: 'auto',
  antialias: true,
  quality : 'medium',
  background: '#EEE'
};

var parent = document.getElementById('viewer');
var viewer = pv.Viewer(parent, options);

function setColorForAtom(go, atom, color){
  var view = go.structure().createEmptyView();
  view.addAtom(atom);
  go.colorBy(pv.color.uniform(color), view);
}

var prevPicked = null;

parent.addEventListener('mousemove', function(event){
  var rect = viewer.boundingClientRect();
  var picked = viewer.pick({ x : event.clientX - rect.left,
        y : event.clientY - rect.top });
  if (prevPicked !== null && picked !== null &&
    picked.target() === prevPicked.atom){
    return;
  }
  if (prevPicked !== null){
    setColorForAtom(prevPicked.node, prevPicked.atom, prevPicked.color);
  }
  if (picked !== null){
    var atom = picked.target();
    var index = atom.residue().num();
    document.getElementById('changes').innerHTML = (index+1) + changes[index]
    var color = [0,0,0,0];
    picked.node().getColorForAtom(atom, color);
    prevPicked = { atom : atom, color : color, node : picked.node() };
    setColorForAtom(picked.node(), atom, 'green');
  }
  else{
    document.getElementById('changes').innerHTML = '&nbsp;';
    prevPicked = null;
  }
  viewer.requestRedraw();
});

function color_substitutions() {
  return new pv.color.ColorOp(function(atom, out, index) {
    annotation = annotations[atom.residue().num()]
    if (annotation == 2) {
      out[index+0] = 1.0; out[index+1] = 0.0;
      out[index+2] = 0.0; out[index+3] = 1.0;
    } else if (annotation == 1) {
      out[index+0] = 0.4; out[index+1] = 0.4;
      out[index+2] = 0.4; out[index+3] = 1.0;
    } else if (annotation == 0) {
      out[index+0] = 1.0; out[index+1] = 1.0;
      out[index+2] = 1.0; out[index+3] = 1.0;
    }
  });
}

function loadSamplePDB() {
  <?php
  	$pdbArray = explode('_', $row['pdbID']);
  	$pdbFile = $pdbArray[0];
  	$pdbChain = $pdbArray[1];
	echo "pv.io.fetchPdb('./pdbs/" . $pdbFile . ".pdb', function(structure) {\n";
	echo "var chain = structure.select({chain: '" . $pdbChain . "'})\n";
  ?>
    viewer.cartoon('protein', chain, { color: color_substitutions() });
    sideChain = chain.select({rindices: indices})
    viewer.ballsAndSticks('sideChain', sideChain);
    viewer.fitTo(chain);
  });
}
window.onresize = function(event) {
  viewer.fitParent();
}
document.addEventListener('DOMContentLoaded', loadSamplePDB);

var rootDiv = document.getElementById("msa");

var opts = {
  el: rootDiv,
  <?php echo 'importURL: "./fasta/' . $famMapID . '.fasta",'; ?>
  vis: {
    labelId: false
  }
};
var m = msa(opts);
m.g.zoomer.set("alignmentHeight", 50) // modifies the default height
</script>
