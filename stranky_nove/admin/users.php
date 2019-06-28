<?php
    include_once('../include/layout_admin/head.php');

    if (isset ($_POST['btn_addUser'])) 
    {

    }

?>  
    <h2>Správa uživatelů</h2>

    <table>
        <tr>
            <th>ID</th>
            <th>Uživ. jméno</th>
            <th>Role</th>
            <th>Úlohy</th>
        </tr>
        <?php

        ?>
    </table>

    <hr>
    <h3>Nový uživatel</h3>
    <form method="post">
        Uživatelské jméno: 
        <input name="userName" type="text">
        <br>
        Heslo: 
        <input name="password" type="password">
        <br>
        Role: 
        <select name="role">
            <option value="admin">Administrátor</option>
            <option value="user">Uživatel</option>
        </select>
        <br>
        <button type="button" id="btn_addUser">Přidat uživatele</button>   
    </form>

<?php 
    include ('../include/layout_admin/foot.php');
?>

