# Schubmult Embed — WordPress plugin

Adds a `[schubmult]` shortcode that embeds the schubmult Flask web widget
(hosted separately) via an `<iframe>` on any WordPress page or post.

## Install

1. ZIP this folder so the archive contains `schubmult-embed/schubmult-embed.php`
   at its root (WordPress's **Upload Plugin** button only accepts `.zip`, not a
   bare `.php` file):

   ```bash
   cd wordpress-plugin
   mkdir -p /tmp/schubmult-embed/schubmult-embed
   cp schubmult-embed.php /tmp/schubmult-embed/schubmult-embed/
   (cd /tmp/schubmult-embed && zip -r /path/to/wordpress-plugin/schubmult-embed.zip schubmult-embed)
   ```

2. In WordPress admin: **Plugins → Add New → Upload Plugin**, choose
   `schubmult-embed.zip`, click **Install Now**, then **Activate**.
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
[schubmult flavor="groth"]
[schubmult flavor="py" height="640" width="100%"]
[schubmult flavor="double"]
[schubmult flavor="groth_q_double"]
```

`flavor` preselects the kernel: `py`, `groth`, `double`, `groth_double`, `q`, `groth_q`,
`q_double`, `groth_q_double`.

## Server-side configuration

The Flask app at `/embed` must allow your WordPress origin to embed it.
On your hosted app set the environment variable:

```
SCHUBMULT_ALLOWED_ORIGINS=https://your-wordpress-site.com,https://www.your-wordpress-site.com
```

Without this, browsers will block the iframe with a CSP `frame-ancestors`
violation.

## Troubleshooting a blank widget

If the page area is blank after reactivating your site/app, check these first:

1. Open your embed URL directly in a browser, e.g.
   `https://YOURUSER.pythonanywhere.com/embed`.
   - If this fails, fix the PythonAnywhere app first.
2. In WordPress, go to **Settings → Schubmult Embed** and re-save the embed URL.
   - Reactivation/migration can clear options.
3. In PythonAnywhere, verify your WSGI file still sets:
   - `SCHUBMULT_ALLOWED_ORIGINS=https://your-site.com,https://www.your-site.com`
   - Then click **Reload** in the Web tab.
4. In browser DevTools (Console), look for
   `Refused to frame ... because an ancestor violates Content Security Policy`.
   - That means your allowed-origins list is missing or incorrect.
