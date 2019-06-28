<?php
  session_start();

  $path=$_SERVER['DOCUMENT_ROOT'] . '/mtsm';
  
  require_once($path.'/include/functions.php');
  require_once($path.'/include/database.php');

  header('Content-Type: text/html; charset=utf-8');
?>

<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>HPC FIT | Administrace</title>
  <link rel="stylesheet" href="../style.css">
</head>

<body>
  <div id="container">
    <div id="containerHead">
  	  <h1>Výzkumná skupina Vysoce náročné výpočty - Administrace</h1>
    </div>
    
    <?php
        include($path.'/include/layout_admin/menu.php');
    ?>

    <div id="containerContent">
      
   		
