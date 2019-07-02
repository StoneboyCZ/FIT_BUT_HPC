<?php
  session_start();

  $path=$_SERVER['DOCUMENT_ROOT'] . dirname($_SERVER['PHP_SELF']);

  require_once($path.'/include/functions.php');
  require_once($path.'/include/database.php');

  header('Content-Type: text/html; charset=utf-8');

  

?>

<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>HPC FIT</title>
  <link rel="stylesheet" href="style.css">
  <script src="jquery-3.4.1.min.js"></script>
</head>

<body>
  <div id="container">
    <div id="containerHead">
  	  <h1>Výzkumná skupina Vysoce náročné výpočty</h1>
    </div>
    
    <?php
        include('include/layout/menu.php');
    ?>

    <div id="containerContent">
      
   		
