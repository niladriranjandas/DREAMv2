<?php
$pathLen = 0;

function prePad($level)
{
  $ss = "";

  for ($ii = 0;  $ii < $level;  $ii++)
  {
    $ss = $ss . "|&nbsp;&nbsp;";
  }

  return $ss;
}

function myScanDir($dir, $level, $rootLen)
{
  global $pathLen;

  if ($handle = opendir($dir)) {

    $allFiles = array();

    while (false !== ($entry = readdir($handle))) {
      if ($entry != "." && $entry != "..") {
        if (is_dir($dir . "/" . $entry))
        {
          $allFiles[] = "D: " . $dir . "/" . $entry;
        }
        else
        {
          $allFiles[] = "F: " . $dir . "/" . $entry;
        }
      }
    }
    closedir($handle);

    natsort($allFiles);

    $preventry = "";
    echo "<table border=1 id=\"table_detail\">";
    echo "<tr>";
    echo "<td>";
    foreach($allFiles as $value)
    {
      $displayName = substr($value, $rootLen + 4);
      $fileName    = substr($value, 3);
      $linkName    = str_replace(" ", "%20", substr($value, $pathLen + 3));
      if (is_dir($fileName)) {        
        //echo prePad($level) . $linkName . "<br>\n";        
        echo "<br>\n";
        echo "<table border=1 id=\"table_detail\">";        
        $tmp = substr($linkName,1);
        echo "<tr onclick=\"showHideRow('" . $tmp  . "');\">";
        echo "<td>";
        echo "<p><center>" . $tmp  . "</center></p>";
        echo "</td>";
        echo "</tr>";
        echo "<tr id=\"" . $tmp  . "\" class=\"hidden_row\">";
        echo "<td>";
        myScanDir($fileName, $level + 1, strlen($fileName));
        echo "</td>";
        echo "</tr>";
        echo "</table>";        
      } else {        
        //echo prePad($level) . "<a href=userpdbs" . $linkName . ">" . $displayName . "</a><br>\n";
        //echo "<tr>";        
        $arr_fields = explode('_',trim($displayName));
        $check_name = $arr_fields[3];
        //echo "$check_name";
        if (strcmp($preventry, $check_name) != 0) {              
              $preventry = $check_name;
              echo "<br>\n";
        }
        $ext =  substr($displayName,-3);
        if (strcmp($ext,"pdb") == 0) {
          echo "PDB coordinate file:" . "&nbsp;&nbsp;<a href=." . $linkName . ">" . $displayName . "</a><br>";
        }
        elseif (strcmp($ext,"pdf") == 0) {
          echo "Report file:" . "&nbsp;<a href=." . $linkName . "><span style=\"padding-left:68px;\"></span>" . $displayName . "</a><br>";
        }
        else {
          echo " ";

        }


        //echo "<td>";
        //$ext =  substr($displayName,-3);
        //if (strcmp($ext,"pdb") == 0) {
        //  echo "<p>PDB coordinate file:</p>";
        //}
        //elseif (strcmp($ext,"pdf") == 0) {
        //  echo "<p>Report file:</p>";
        //}
        //else {
        //  echo "<p>Extra file:</p>";

        //}
        //echo "</td>";
        //echo "<td>";
        //echo "<p>" . "<a href=userpdbs" . $linkName . ">" . $displayName . "</a><br>"  . "</p>";
        //echo "</td>";

        //echo "</tr>";
      }
    }
    echo "</td>";
    echo "</tr>";
    echo "</table>";
  }
}

?><!DOCTYPE HTML>
<html>
  
<head>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.2/js/bootstrap.min.js"></script>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.2/css/bootstrap.min.css">
    <link rel="stylesheet" type="text/css" href="https://use.fontawesome.com/releases/v5.6.3/css/all.css">
  
    <script type="text/javascript">
        function showHideRow(row) {
            $("#" + row).toggle();
        }
    </script>
  
    <style>
        body {
            margin: 0 auto;
            padding: 0px;
            text-align: center;
            width: 100%;
            font-family: "Myriad Pro", 
                "Helvetica Neue", Helvetica, 
                Arial, Sans-Serif;
        }
  
        #wrapper {
            margin: 0 auto;
            padding: 0px;
            text-align: center;
            width: 995px;
        }
  
        #wrapper h1 {
            margin-top: 50px;
            font-size: 45px;
            color: #585858;
        }
  
        #wrapper h1 p {
            font-size: 20px;
        }
  
        #table_detail {
            width: 800px;
            text-align: left;
            border-collapse: collapse;
            color: #2E2E2E;
            border: #A4A4A4;
        }
  
        #table_detail tr:hover {
            background-color: #F2F2F2;
        }
  
        #table_detail .hidden_row {
            display: none;
        }
    </style>
</head>


<body>
<h3>DREAM (Distance Restraint and Energy Assisted Modeling)</h3>
<div id="wrapper">
<?php
  $root = '.';

  $pathLen = strlen($root);  
  myScanDir($root, 0, strlen($root)); ?>
</div>
</body>

</html>
