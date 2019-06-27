<?php
    include_once("global.php");

    if (!isset($_SESSION['cookie'])) // pokud uzivatel neni prihlasen
    {
        // bude presmerovan na hlavni stranku
        header ("Location: http://" . $_SERVER['HTTP_HOST'] . dirname($_SERVER['PHP_SELF']) . "/main.php");
        ob_end_clean();
        exit();
    }
    else // jinak odhlasime uzivatele
    {
        $_SESSION = array();  // smazeme promenne
        session_destroy(); // ukoncime session
        setcookie(session_name(), '', time()-300, '/', '', 0); // smazeme cookie

    }


?>

<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Excel@FIT IS</title>
</head>

<body>
    <h1>Cookie smazána.</h1>

    <a href="index.php">Zpět na zadání cookie</a> 
</body>

</html> 