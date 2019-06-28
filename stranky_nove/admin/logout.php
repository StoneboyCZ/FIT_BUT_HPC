<?php
    #include_once("global.php");

    if (!isset($_SESSION['cookie'])) // pokud uzivatel neni prihlasen
    {
        // bude presmerovan na hlavni stranku
        header ("Location: http://" . $_SERVER['HTTP_HOST'] . dirname($_SERVER['PHP_SELF']) . "/index.php");
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
    <title>FIT HPC</title>
</head>

<body>
    <h2>Odhlášení úspěšné.</h2>

    <a href="../index.php">Zpět na stránky skupiny</a>
    <a href="index.php">Znovu přihlásit do administrace</a>  
</body>

</html> 