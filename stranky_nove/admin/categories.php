<?php
    include_once('../include/layout_admin/head.php');

    if (isset ($_POST['btn_addCategory'])) 
    {
        $id=$categories->count()+1; // new user ID (number of user +1)
        
        // create new item
        $item = $categories->get($id);

        // set values
        $item->id = $id;
        $item->name = $_POST['categoryName'];
        $item->order = $_POST['order'];
        $item->status = 0;
        
        $item->save();

        // block resubmission of POST data
        header("Location: ".$_SERVER['PHP_SELF']);
    }
    
    if (isset ($_POST['btn_removeCategory']))
    {
        if ($categories->has($_POST['categoryId']))
        {
            $item = $categories->get($_POST['categoryId']);
            
            // the user is no longer valid
            $item->status = 1; 
            $item->save();       
        }
        
        // block resubmission of POST data
        header("Location: ".$_SERVER['PHP_SELF']);    
    }



?>  
    <h2>Správa kategorií</h2>
    <table>
        <tr>
            <th>ID</th>
            <th>Název</th>
            <th>Pořadí</th>
            <th>Úlohy</th>
        </tr>
        <?php
            $validCategories = $categories->where('status','=',0)->resultDocuments();
            foreach($validCategories as $category) 
            {
                //echo $category;
                if ($category->status == 0)
                {
                    echo '<tr>';
                    echo '<td>'.$category->id.'</td>';
                    echo '<td>'.$category->name.'</td>';
                    echo '<td>'.$category->order.'</td>';
                    echo '<td><form method="post" action="'.$_SERVER['PHP_SELF'].'">
                    <input type="hidden" id="categoryId" name="categoryId" value="'.$category->id.'">
                    <input type="submit" name="btn_updateCategory" value="Upravit kategorii">
                    <input type="submit" name="btn_removeCategory" value="Odstranit kategorii"></td>
                    </form>'; 
                    //echo '<td>Editovat    '.'<a href="#" onclick="removeUser('.$user->id.')" id="remUser_'.$user->id.'">Odebrat</a></td>';
                    echo '</tr>';
                }
                    
            }
        ?>
    </table>

    <hr>
    <h3>Nová kategorie</h3>
    <form method="post" action="<?php echo $_SERVER['PHP_SELF']; ?>">
        Název kategorie: 
        <input name="categoryName" type="text">
        <br>
        Pořadí: 

        <select name="order">
            <?php
                for($i=0; $i<count($validCategories)+1;$i++)
                {
                    echo '<option value="'.$i.'">'.$i.'</option>';
                }    
            ?>
        </select>
        <br>
        <input type="submit" value="Přidat kategorii" name="btn_addCategory"> 
    </form>



<?php 
    include ('../include/layout_admin/foot.php');
?>

