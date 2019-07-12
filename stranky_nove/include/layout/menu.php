    <div id="containerMenu">
        <ul>
          <li><a href="index.php">Úvodní stránka</a></li>
        <?php
          $validCategories = $categories->where('status','=',0)->resultDocuments();
          foreach($validCategories as $category)
          {
            echo '<li class="menu_label">'.$category->name.'</li>';  
          } 
        ?>
        </ul>
      <ul>
      <li><a href="admin/index.php">Administrace</a></li>
        </ul> 

    </div>