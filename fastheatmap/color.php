<html>
<body>

<?php
// define variables and set to empty values
$colorErr = $lowColorErr   = "";
$color = $lowColor    = "";

if ($_SERVER["REQUEST_METHOD"] == "POST") {
  if (empty($_POST["highcolor"])) {
    $colorErr = "Name is required";
  } else {
    $color = test_input($_POST["highcolor"]);
    // check if name only contains letters and whitespace
    if (!preg_match("/^[a-zA-Z ]*$/",$color)) {
      $colorErr = "Only letters and white space allowed"; 
    }
  }
  
  if (empty($_POST["lowcolor"])) {
    $lowColorErr = "Email is required";
  } else {
    $lowColor = test_input($_POST["lowcolor"]);
    // check if e-mail address is well-formed
   
  }
 

 

  
}

function test_input($data) {
  $data = trim($data);
  $data = stripslashes($data);
  $data = htmlspecialchars($data);
  return $data;
}
?>
<?php
return $color;
?>
</body>
</html>
