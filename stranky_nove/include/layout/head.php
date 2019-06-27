<?php
  session_start();

  require_once('include/functions.php');
  require_once('include/database.php');

  header('Content-Type: text/html; charset=utf-8');
?>

<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>HPC FIT</title>
  <link rel="stylesheet" href="style.css">
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
      
   		
