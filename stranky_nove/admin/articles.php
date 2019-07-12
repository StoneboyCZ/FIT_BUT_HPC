<?php
    include_once('../include/layout_admin/head.php');

    if (isset ($_POST['btn_addArticle'])) 
    {
        $id=$articles->count()+1; // new user ID (number of user +1)
        
        // create new item
        $item = $articles->get($id);

        // set values
        $item->id = $id;
        $item->name = $_POST['articleName'];
        $item->category = $_POST['categoryName'];
        
        // get the category for the item
        $catItem = $categories->get($item->category);

        // add the new article to the end of the category
        $catCount = count($catItem->articles);
        $item->order = $catCount;

        // insert article into the category
        array_push($catItem->articles,$id);         
        $catItem->save();

        $item->status = 0;
        
        $item->save();

        // block resubmission of POST data
        header("Location: ".$_SERVER['PHP_SELF']);
    }
    
    if (isset ($_POST['btn_removeArticle']))
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
    <h2>Správa článků</h2>
    <table>
        <tr>
            <th>ID</th>
            <th>Název</th>
            <th>Kategorie</th>
            <th>Pořadí</th>
            <th>Úlohy</th>
        </tr>
        <?php
            $validArticles = $articles->where('status','=',0)->resultDocuments();
            foreach($validArticles as $i=>$article) 
            {
                echo '<tr>';
                echo '<td>'.$article->id.'</td>';
                echo '<td>'.$article->name.'</td>';
                echo '<td>'.$article->category.'</td>';
                echo '<td>'.$article->order.'</td>';
                echo '<td><form method="post" action="'.$_SERVER['PHP_SELF'].'">
                <input type="hidden" id="articleId" name="categoryId" value="'.$article->id.'">
                <input type="submit" name="btn_updateArticle" value="Upravit článek">
                <input type="submit" name="btn_removeArticle" value="Odstranit článek"></td>
                </form>'; 
                echo '</tr>';
            }
        ?>
    </table>

    <hr>
    <h3>Nový článek</h3>
    <form method="post" action="<?php echo $_SERVER['PHP_SELF']; ?>">
        Název článku: 
        <input name="articleName" type="text">
        <br>
        Kategorie 
        <select name="categoryName">
            <?php
                $validCategories = $categories->where('status','=',0)->resultDocuments();
                foreach($validCategories as $category)
                {
                    echo '<option value="'.$category->id.'">'.$category->name.'</option>';
                }
            ?>
        </select>
        <br>
        <input type="submit" value="Přidat článek" name="btn_addArticle"> 
    </form>



<?php 
    include ('../include/layout_admin/foot.php');
?>

