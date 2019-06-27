<?php 

include ('layout/head.php');

if (isset($_SESSION['cookie'])) // pokud jsme prihlaseni 
{?> 
    <div id="containerContent">
    <h2>Vyberte sekci z menu.</h2>
    </div>   
    
<?php }
//else include ('error.php'); 

include ('layout/foot.php');

?>
