# Schubmult Embed — WordPress plugin

Adds a `[schubmult]` shortcode that embeds the schubmult Flask web widget
(hosted separately) via an `<iframe>` on any WordPress page or post.

## Install

1. ZIP this folder so the archive contains `schubmult-embed/schubmult-embed.php`
   at its root, or just upload the `.php` file directly.
2. In WordPress admin: **Plugins → Add New → Upload Plugin**, choose the file,
   click **Install Now**, then **Activate**.
3. Go to **Settings → Schubmult Embed** and paste your hosted app URL,
   e.g. `https://yourname.pythonanywhere.com/embed`.

## Use

In any post or page (Gutenberg or Classic editor) insert:

```
[schubmult]
```

Optional attributes:

```
[schubmult height="500"]
[schubmult flavor="py" height="640" width="100%"]
[schubmult flavor="double"]
```

## Server-side configuration

The Flask app at `/embed` must allow your WordPress origin to embed it.
On your hosted app set the environment variable:

```
SCHUBMULT_ALLOWED_ORIGINS=https://your-wordpress-site.com,https://www.your-wordpress-site.com
```

Without this, browsers will block the iframe with a CSP `frame-ancestors`
violation.
