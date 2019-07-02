<?php
    include_once('../include/layout_admin/head.php');

    if (isset ($_POST['btn_addUser'])) 
    {
        $id=$users->count()+1; // new user ID (number of user +1)
        
        // create new item
        $item = $users->get($id);

        // set values
        $item->id = $id;
        $item->userName = $_POST['userName'];
        $item->password = hash('sha256',$_POST['password']);
        $item->role = $_POST['role'];
        $item->status = 0;
        
        $item->save();

        // block resubmission of POST data
        header("Location: ".$_SERVER['PHP_SELF']);
    }

    if (isset ($_POST['btn_removeUser']))
    {
        if ($users->has($_POST['userId']))
        {
            $item = $users->get($_POST['userId']);
            
            // the user is no longer valid
            $item->status = 1; 
            $item->save();       
        }
        
        // block resubmission of POST data
        header("Location: ".$_SERVER['PHP_SELF']);    
    }



?>  
    <script type="text/javascript">
        function removeUser(id)
        {
            
        }
    </script>
    <h2>Správa uživatelů</h2>

    

    <table>
        <tr>
            <th>ID</th>
            <th>Uživ. jméno</th>
            <th>Role</th>
            <th>Úlohy</th>
        </tr>
        <?php
            $results = $users->findAll();
            
            foreach($results as $user) 
            {
                if ($user->status == 0)
                {
                    echo '<tr>';
                    echo '<td>'.$user->id.'</td>';
                    echo '<td>'.$user->userName.'</td>';
                    echo '<td>'.$user->role.'</td>';
                    echo '<td><form method="post" action="'.$_SERVER['PHP_SELF'].'">
                    <input type="hidden" id="userId" name="userId" value="'.$user->id.'">
                    <input type="submit" name="btn_removeUser" value="Odstranit uživatele">
                    </form>'; 
                    //echo '<td>Editovat    '.'<a href="#" onclick="removeUser('.$user->id.')" id="remUser_'.$user->id.'">Odebrat</a></td>';
                    echo '</tr>';
                }
                    
            }
        ?>
    </table>

    <hr>
    <h3>Nový uživatel</h3>
    <form method="post" action="<?php echo $_SERVER['PHP_SELF']; ?>">
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
        <input type="submit" value="Přidat uživatele" name="btn_addUser"> 
    </form>



<?php 
    include ('../include/layout_admin/foot.php');
?>

